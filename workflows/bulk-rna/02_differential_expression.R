# Script: 02_differential_expression.R
# Purpose: Differential expression analysis with DESeq2
# Workflow: bulk-rna
# Required inputs: Bulk_batch_count.csv; group.txt
# Main outputs: Bulk_batch_differ.csv
# Prerequisites: Requires batch-corrected counts from step 01.

rm(list = ls())

library(DESeq2)
expr_matrix <- read.csv("Bulk_batch_count.csv", row.names = 1)
expr_matrix <- round(expr_matrix)

metadata <- read.delim('./group.txt',header = T, stringsAsFactors = F)
rownames(metadata) <- metadata$sample
metadata <- metadata[colnames(expr_matrix),]
metadata$group <- factor(metadata$group, levels = c("Low", "High"))

colnames(expr_matrix)
condition <- factor(metadata$group, levels = c("Low","High"))    
coldata <- data.frame(row.names = colnames(expr_matrix), condition)
coldata   

dds <- DESeqDataSetFromMatrix(countData=expr_matrix, colData=coldata, design=~condition)

keep <- rowSums(counts(dds) >= 10) >= 3  
dds <- dds[keep, ]  
dds1 <- DESeq(dds)   
resultsNames(dds1)   
dds1$condition     
res <- results(dds1, pAdjustMethod = "BH")
summary(res)       

table(res$padj < 0.05)       
res <- res[order(res$padj),]  
res$Gene <- rownames(res)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)

write.csv(resdata,file = "Bulk_batch_differ.csv")   