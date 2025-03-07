# RNA Seq workflow: Differential gene expression analysis
## Differential gene expression analysis based on the negative binomial distribution using Bioconductor - DESeq2 R package
### 1. Installing and loading required Bioconductor R packages
Install BiocManager, if required
```{r}
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install(version = "3.19")
}
```
Install and load all bioconductor packages
```{r}
packages <- c("DESeq2",
              "pheatmap",
              "RColorBrewer",
              "tidyverse",
              "apeglm")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]

if(length(new.packages))
  BiocManager::install(new.packages)

invisible(lapply(packages, library, character.only=TRUE))
```
### 2. Setting working directory and loading the data
Set the directory having the current script as 'working directory' and create an output file.
```{r}
if (!require("rstudioapi")) {
  install.packages("rstudioapi")
  library(rstudioapi)
}
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ifelse(!dir.exists(file.path("./Output")),
  dir.create(file.path("./Output")),
  "Directory Exists"
)
list.files()
```
Load the data file
```{r}
rawcount <- read.table("./Input/rawcounts.csv",
  sep = ",",
  header = TRUE,
  row.names = 1
)

head(rawcount) # Check column names
dim(rawcount) # Check row count

```
### 3. Data Preparation
Define and create metadata table
```{r}
# Define metadata
sample_type <- c(
  "smoc2_oe",
  "smoc2_oe",
  "smoc2_oe",
  "smoc2_oe",
  "smoc2_oe",
  "smoc2_oe",
  "smoc2_oe"
)
condition <- c(
  "fibrosis",
  "fibrosis",
  "fibrosis",
  "fibrosis",
  "normal",
  "normal",
  "normal"
)
row_names <- c(
  "smoc2_fibrosis1",
  "smoc2_fibrosis2",
  "smoc2_fibrosis3",
  "smoc2_fibrosis4",
  "smoc2_normal1",
  "smoc2_normal3",
  "smoc2_normal4"
)
color <- c(
  "salmon1",
  "salmon1",
  "salmon1",
  "salmon1",
  "cadetblue1",
  "cadetblue1",
  "cadetblue1"
)

# Create metadata datadrame
metadata <- data.frame(
  row.names = TRUE,
  row_names,
  sample_type,
  condition,
  color
)

metadata # View the table
```

Reorder the columns of rawcount table to match the row order of metadata
```{r}
# Check the data order with respect to metadata
all(colnames(rawcount) == rownames(metadata))

# Define a function to reorder data columns to match metadata
fn <- match(
  rownames(metadata),
  colnames(rawcount)
)

# Reorder the data columns
reordered_rawcounts <- rawcount[, fn]

# Confirm the data order (all() function can also be used)
identical(
  colnames(reordered_rawcounts),
  rownames(metadata)
)
```

Visualize rawcount data
```{r}
barplot(colSums(reordered_rawcounts) / 1000000,
  col = metadata$color,
  ylab = "Million counts",
  cex.names = 1,
  las = 2,
  ylim = c(0, 30)
)
```

Normalize gene counts
```{r}
# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = reordered_rawcounts,
  colData = metadata,
  design = ~condition
)

# Determine the size factors
dds <- estimateSizeFactors(dds)

# Extract the normalized counts
normalized_count <- counts(dds,
  normalized = TRUE
)
```

Visualize normalized counts
```{r}
# Log transformation of the the normalized counts
vsd <- vst(dds, blind = TRUE)

# Plot the PCA of PC1 and PC2
plotPCA(vsd,
  intgroup = "condition"
)

vsd_mat <- assay(vsd) # Extract the transformed counts matrix
vsd_cor <- cor(vsd_mat) # Compute the correlation between samples

# Plot the heatmap for corelation values
pheatmap(vsd_cor)
```

Estimate size factors and dispersion followed by fitting negative binomial model (Wald test)
```{r}
dds <- DESeq(dds)
```

Check default order of comparison selected, which is indicated by "condition_test_vs_reference". If the wrong group is set as reference, re-level and run wald test again to change the order.
```{r}
# Check for comparison order; expectation: condition_fibrosis_vs_normal
resultsNames(dds)

# relevel to set normal as reference,
dds$condition <- relevel(dds$condition, ref = "normal")
# re-run waldTest
dds <- nbinomWaldTest(dds)
# check comparison order again
resultsNames(dds)
```

Explore the fit of the data to the negative binomial model by plotting the dispersion estimates
```{r}
plotDispEsts(dds)
```

Extract the log2 fold changes of the fibrosis group relative to normal
Setting a threshold to fold change (not always prefered) will return lesser but biologically relevant DEGs 
```{r}
# Extract results
res <- results(dds,
  alpha = 0.05,
  lfcThreshold = 0.32
) # 0.32 lfcThreshold for 1.25 fold change

# MA plot before shrinking
plotMA(res,
  ylim = c(-8, 8)
)
```

Improve the estimated log2 fold change by shrinking the results
Note: This might throw a warning if results have transcripts with baseMean=0.
```{r}
# Shrink the log2 fold changes
res <- lfcShrink(dds,
  type = "apeglm",
  coef = 2,
  res = res
)
# MA plot after shrinking
plotMA(res,
  ylim = c(-8, 8)
)

```

Extract and save the significant DEGs from the results with padj < 0.05
```{r}
summary(res) # Overview of the results

# Generate result table
deg <- data.frame(res)
# Create subset of results having only the significant DEGs
deg_sig <- subset(
  deg,
  padj < 0.05
)

# Save results
write.table(deg_sig,
  "Output/smoc2_DEG_sig.txt",
  sep = "\t",
  row.names = TRUE
)
```

Visualize the results:

Volcano plot of significant DEGs, where red and green dots represent down-regulated and up-regulated genes, respectively, with the fold change of greater than 1.
```{r}
de <- deg[complete.cases(deg), ] # Filter out NA
# Add new column to define differential expression
de$diffexpressed <- "NO"
# Define up-regulated and down-regulated genes
de$diffexpressed[de$log2FoldChange > 1 & de$padj < 0.05] <- "UP"
de$diffexpressed[de$log2FoldChange < -1 & de$padj < 0.05] <- "DOWN"
# Asign colors
mycolors <- c("red", "green", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

# Volcano plot
ggplot(
  data = de,
  aes(
    x = log2FoldChange,
    y = -log10(pvalue),
    col = diffexpressed
  )
) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values = mycolors) +
  geom_vline(
    xintercept = c(-1, 1),
    col = "black"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    col = "black"
  )
```

Heatmap of normalized counts of significant DEGs
```{r}
# Generate a table with normalized counts
normalized_count_table <- data.frame(normalized_count)

# Subset the table having only significant DEGs
sig_norm_counts <- normalized_count_table[rownames(deg_sig), ]

# Plot heatmap
pheatmap(sig_norm_counts,
  cluster_rows = TRUE,
  show_rownames = FALSE,
  scale = "row"
)
```
