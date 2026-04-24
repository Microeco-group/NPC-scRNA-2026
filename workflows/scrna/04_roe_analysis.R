# Script: 04_roe_analysis.R
# Purpose: Relative enrichment (Ro/e) analysis across groups and samples
# Workflow: scRNA
# Required inputs: Annotated Seurat object metadata with celltype/group/sample information.
# Main outputs: Roe_result.rds
# Prerequisites: Requires a curated Seurat object with final annotations.

rm(list = ls())
library(tictoc)
library(tidyverse)
library(sscVis)
library(Seurat)
library(tidyverse)
library(readr)
library(qs)
library(BiocParallel)
library(data.table)
library(plyr)

Idents(scRNA) <- "celltype"
data <- scRNA@meta.data
colnames(data)

Roe <- calTissueDist(data,
                     byPatient = F,
                     colname.cluster = "celltype", 
                     colname.patient = "sample",
                     colname.tissue = "group",
                     method = "chisq", 
                     min.rowSum = 0) 
Roe
saveRDS(Roe, "Roe_result.rds")