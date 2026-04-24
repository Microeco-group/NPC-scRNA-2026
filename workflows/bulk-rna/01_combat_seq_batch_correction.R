# Script: 01_combat_seq_batch_correction.R
# Purpose: ComBat-seq batch correction for bulk RNA-seq counts
# Workflow: bulk-rna
# Required inputs: norm_data_count.csv; group.txt with sample, group, and dataset columns.
# Main outputs: Bulk_batch_count.csv
# Prerequisites: Run from this workflow directory so relative paths resolve correctly.

rm(list = ls())
library(DESeq2)
library(sva)

metadata <- read.delim('group.txt',header = T, stringsAsFactors = F)
rownames(metadata) <- metadata$sample
metadata <- metadata[colnames(norm_data),]
metadata$group <- factor(metadata$group, levels = c("Low", "High"))

norm_data <- read.csv("norm_data_count.csv")
colnames(norm_data)
condition <- factor(metadata$group, levels = c("Low","High"))    
coldata <- data.frame(row.names = colnames(norm_data), condition)
coldata  

dds <- DESeqDataSetFromMatrix(countData=norm_data, colData=coldata, design=~condition)
keep <- rowSums(counts(dds) >= 10) >= 3  
dds <- dds[keep, ]    
dds1 <- DESeq(dds)   
normalized_counts <- counts(dds1, normalized = TRUE)

#correct
norm_data <- as.matrix(normalized_counts)
batch <- metadata$dataset
mod <- metadata$group
combat_data <- ComBat_seq(norm_data, batch = batch)

min(combat_data)
max(combat_data)

write.csv(combat_data, "Bulk_batch_count.csv")