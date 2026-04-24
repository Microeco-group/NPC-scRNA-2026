# Script: 12_nmf_program_discovery.R
# Purpose: NMF program discovery in epithelial cells
# Workflow: scRNA
# Required inputs: EpithelialCells_clean.rds
# Main outputs: res_nmf_sampleX.rds
# Prerequisites: Requires a malignant epithelial-cell subset.

rm(list=ls())
library(Seurat)
library(NMF)      
library(dplyr)    
library(tidyr)    

EpithelialCells_clean <- readRDS("EpithelialCells_clean.rds")
Target_seurat <- subset(EpithelialCells_clean, subset = celltype2 == "MalignantCells")

DefaultAssay(Target_seurat) <- "RNA"
table(Target_seurat@meta.data$sample)

rank_k <- 2:10
exp_matrix <- Target_seurat@assays$RNA@scale.data
res_nmf <- nmf(exp_matrix, rank_k, seed = 10, method = "snmf/r", .options = "vtp30")
saveRDS(res_nmf,"res_nmf_sampleX.rds")
