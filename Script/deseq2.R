#-----Libraries-----
#Install BiocManager as per Bioconductor instructions
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")

#Install required packages
BiocManager::install(c("DESeq2", "affy","pheatmap","RColorBrewer","tidyverse", "apeglm"))

#Load libraries
library(DESeq2)
library(affy)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(apeglm)

#-----Importing, defining and cleaning the data-----
#Set working directory
setwd("C:/Users/anush/Documents/Academic/DESeq2_fibrosis_smoc2/")
list.files()

#Explore rawcount table
rawcount<-read.table("Input/fibrosis_smoc2_rawcounts.csv",sep = ",", header = T, row.names = 1)
head(rawcount) #Check column names
dim(rawcount) #Check row count

#Define metadata 
sample_type <- c('smoc2_oe','smoc2_oe','smoc2_oe','smoc2_oe','smoc2_oe','smoc2_oe','smoc2_oe')
condition <- c("fibrosis","fibrosis","fibrosis","fibrosis","normal","normal","normal")
row_names <- c('smoc2_fibrosis1','smoc2_fibrosis2','smoc2_fibrosis3','smoc2_fibrosis4','smoc2_normal1','smoc2_normal3','smoc2_normal4')
defcolor <- c("fibrosis"="salmon1", "normal"="cadetblue1")

metadata <- data.frame(row.names = TRUE, row_names, sample_type,condition)
metadata$color <- defcolor[as.vector(metadata$condition)]

head(metadata) #Check table structure
table(metadata$condition) #Count number of samples per condition

# Reorder the columns of rawcount table to match metadata rows
reorder_fn <- match(rownames(metadata),colnames(rawcount))
reordered_rawcounts <- rawcount[,reorder_fn]
if(identical(colnames(reordered_rawcounts), rownames(metadata))==FALSE) 
  {print("sample order of the dataset do not match with rows of metadata")}

#-----visualize rawcount data---------
png(filename = "Reads_per_sample.png")
barplot(colSums(reordered_rawcounts)/1000000, 
        main="Total number of reads per sample", 
        col = metadata$color, 
        ylab = "Million counts", 
        cex.names = 1, 
        las = 2, 
        ylim = c(0,30))
dev.off()

#----------Differential gene expression analysis----------
# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = reordered_rawcounts,
                                  colData =  metadata,
                                  design = ~ condition)

#Determine the size factors to use for normalization
dds <- estimateSizeFactors(dds)

# Extract the normalized counts
normalized_count <- counts(dds, normalized=TRUE)

#---------Plot normalized data---------
vsd <- vst(dds, blind = TRUE) # Transform the normalized counts 
# Plot the PCA of PC1 and PC2
png(filename = "PCA_normalized_counts.png")
plotPCA(vsd, intgroup='condition')
dev.off()
vsd_mat <- assay(vsd) # Extract the matrix of transformed counts
vsd_cor <- cor(vsd_mat) # Compute the correlation values between samples
# Plot the heatmap
png(filename = "Heatmap_normalized_counts.png")
pheatmap(vsd_cor,
         main = "Heatmap of normalized counts")
dev.off()

#----Dispersion----
# Plot dispersions
dds_dis<-DESeq(dds)
png(filename = "Dispersion_counts.png")
plotDispEsts(dds_dis)
dev.off()


#--------results---------
# Extract results
res <- results(dds_dis,
                     #contrast = c("condition", "fibrosis", "normal"), 
                     alpha = 0.05, 
                     lfcThreshold = 0.32) #lfcThreshold = 0.32 corresponds to 1.25 fold change

# Shrink the log2 fold changes
res <- lfcShrink(dds_dis, 
                       #contrast = c("condition", "fibrosis", "normal"), 
                     type = 'apeglm', coef = 2, res = res)

#-------MA plot ---------
# Create MA plot
png(filename = "MAplot.png")
plotMA(res, main='MA-plot', alpha=0.1)
abline(h=c(-1:1), col="red")
dev.off()

#----------filter and save results-----------
summary(res) # Get an overview of the results          
deg <- data.frame(res)# Save results as a data frame
deg_sig <- subset(deg, padj < 0.05)# Subset the results to only return the significant genes

#-----------volcano plot-----------
de <- deg[complete.cases(deg), ]
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > 1 & de$padj < 0.05] <- "UP"
de$diffexpressed[de$log2FoldChange < -1 & de$padj < 0.05] <- "DOWN"
mycolors <- c("red", "green", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  scale_colour_manual(values = mycolors) +
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")

#----heatmap------
normalized_count_table<-data.frame(normalized_count)

# Subset normalized counts to significant genes
sig_norm_counts <- normalized_count_table[rownames(deg_sig), ]

# Plot heatmap
pheatmap(sig_norm_counts, 
         cluster_rows = T, 
         show_rownames = F, 
         scale = 'row')
