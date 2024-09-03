#-----Libraries-----
#Install BiocManager as per Bioconductor instructions
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")

#Install required packages
BiocManager::install(c("DESeq2",
                       "affy",
                       "pheatmap",
                       "RColorBrewer",
                       "tidyverse",
                       "apeglm"))

#Load required packages
library(DESeq2)
library(affy)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(apeglm)

#-----Importing, defining and reordering the data-----
#Set working directory
setwd("C:/workspace/github/differential-expression-analysis/Output")
list.files()

#Explore rawcount table
rawcount <- read.table("C:/workspace/github/differential-expression-analysis/Input/rawcounts.csv",
                       sep = ",",
                       header = TRUE,
                       row.names = 1)

head(rawcount) #Check column names
dim(rawcount) #Check row count

#Define metadata
sample_type <- c("smoc2_oe",
                 "smoc2_oe",
                 "smoc2_oe",
                 "smoc2_oe",
                 "smoc2_oe",
                 "smoc2_oe",
                 "smoc2_oe")
condition <- c("fibrosis",
               "fibrosis",
               "fibrosis",
               "fibrosis",
               "normal",
               "normal",
               "normal")
row_names <- c("smoc2_fibrosis1",
               "smoc2_fibrosis2",
               "smoc2_fibrosis3",
               "smoc2_fibrosis4",
               "smoc2_normal1",
               "smoc2_normal3",
               "smoc2_normal4")
defcolor <- c("fibrosis" <- "salmon1",
              "normal" <- "cadetblue1")

#Create metadata datadrame
metadata <- data.frame(row.names = TRUE,
                       row_names,
                       sample_type,
                       condition)
#Asign colors to conditions
metadata$color <- defcolor[as.vector(metadata$condition)]

head(metadata) #Check table structure
table(metadata$condition) #Count number of samples per condition

# Define reorder parameters to match metadata
fn <- match(rownames(metadata),
            colnames(rawcount))

#Reorder the colums to match metadata
reordered_rawcounts <- rawcount[, fn]

# Check the order
if (identical(colnames(reordered_rawcounts),
              rownames(metadata)) == FALSE)
  print("sample order of the dataset do not match with rows of metadata")

#-----Visualize rawcount data---------
png(filename = "Rawreads_per_Sample.png")
barplot(colSums(reordered_rawcounts) / 1000000,
        col = metadata$color,
        ylab = "Million counts",
        cex.names = 1,
        las = 2,
        ylim = c(0, 30))
dev.off()

#----------Differential gene expression analysis----------
# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = reordered_rawcounts,
                              colData =  metadata,
                              design = ~ condition)

#Determine the size factors
dds <- estimateSizeFactors(dds)

# Extract the normalized counts
normalized_count <- counts(dds,
                           normalized = TRUE)

#---------Plot normalized data---------
vsd <- vst(dds, blind = TRUE) # Transform the normalized counts
# Plot the PCA of PC1 and PC2
png(filename = "PCA_Normalized_Counts.png")
plotPCA(vsd,
        intgroup = "condition")
dev.off()

vsd_mat <- assay(vsd) # Extract the matrix of transformed counts
vsd_cor <- cor(vsd_mat) # Compute the correlation values between samples
# Plot the heatmap
png(filename = "Heatmap_normalized_counts.png")
pheatmap(vsd_cor,
         main = "Heatmap of normalized counts")
dev.off()

# Plot dispersions
dds_dis <- DESeq(dds)
png(filename = "Dispersion_Normalized_Counts.png")
plotDispEsts(dds_dis)
dev.off()

#--------Extract and filter results---------
# Extract results
res <- results(dds_dis,
               alpha = 0.05,
               lfcThreshold = 0.32) #0.32 lfcThreshold for 1.25 fold change

# Shrink the log2 fold changes
res <- lfcShrink(dds_dis,
                 type = "apeglm",
                 coef = 2,
                 res = res)

summary(res) #Overview of the results

#MA plot
png(filename = "MA_plot_Log2FC.png")
plotMA(res,
       alpha = 0.1)
abline(h = c(-1:1),
       col = "red")
dev.off()

#Generate result table
deg <- data.frame(res)
# Create subset of results having only the significant DEGs
deg_sig <- subset(deg,
                  padj < 0.05)

#---------Volcano plot------------
de <- deg[complete.cases(deg), ] #Filter out NA
#Add new column to define differential expression
de$diffexpressed <- "NO"
#Define up-regulated and down-regulated genes
de$diffexpressed[de$log2FoldChange > 1 & de$padj < 0.05] <- "UP"
de$diffexpressed[de$log2FoldChange < -1 & de$padj < 0.05] <- "DOWN"
#Asign colors
mycolors <- c("red", "green", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

#Volcano plot
png(filename = "Volcano_Plot_Log2FC.png")
ggplot(data = de,
       aes(x = log2FoldChange,
           y = -log10(pvalue),
           col = diffexpressed)) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values = mycolors) +
  geom_vline(xintercept = c(-1, 1),
             col = "black") +
  geom_hline(yintercept = -log10(0.05),
             col = "black")
dev.off()

#----Heatmap------
#Generate a table with normalized counts
normalized_count_table <- data.frame(normalized_count)

# Subset the table having only significant DEGs
sig_norm_counts <- normalized_count_table[rownames(deg_sig), ]

# Plot heatmap
png(filename = "DEG_Heatmap.png")
pheatmap(sig_norm_counts,
         cluster_rows = TRUE,
         show_rownames = FALSE,
         scale = "row")
dev.off()

#Save results
write.table(deg_sig,
            "DESeq2_DEG_sig.txt",
            sep = "\t",
            row.names = TRUE)

quit()