# Script: 01_seurat_integration.R
# Purpose: Seurat normalization, integration, and clustering
# Workflow: scRNA
# Required inputs: scRNA.rds; genes_black object available in the session or loaded before execution.
# Main outputs: scRNA_end.rds
# Prerequisites: Run from this workflow directory so relative paths resolve correctly.

rm(list=ls())
library(Seurat)
library(dplyr)
library(tidyverse)
library(SeuratObject)
library(harmony)
library(clustree)

scRNA <- readRDS("scRNA.rds")
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst")
sum(scRNA@assays$RNA@var.features %in% genes_black)
length(scRNA@assays$RNA@var.features)
scRNA@assays$RNA@var.features = scRNA@assays$RNA@var.features[!(scRNA@assays$RNA@var.features %in% genes_black)]
length(scRNA@assays$RNA@var.features)
scRNA <- ScaleData(scRNA,vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mito","percent.rb","S.Score","G2M.Score"))
scRNA <- RunPCA(scRNA, verbose = FALSE)
ElbowPlot(scRNA,ndims = 50)
scRNA <- scRNA %>% RunHarmony("orig.ident",plot_convergence =T)
scRNA <- RunUMAP(scRNA, reduction="harmony",dims=1:30) 
scRNA <- FindNeighbors(scRNA, reduction="harmony", dims = 1:30, verbose = FALSE)   #dims和RunUMAP一致
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 , 1.2, 1.5)) {
  scRNA <- FindClusters(scRNA, 
                              resolution = res, verbose = FALSE)
}
clustree(scRNA)
saveRDS(scRNA,"scRNA_end.rds")