# Script: 03_marker_heatmap.R
# Purpose: Marker-expression heatmap generation
# Workflow: scRNA
# Required inputs: scRNA.rds with curated cell-type markers and annotations.
# Main outputs: Heatmap figures written by the script.
# Prerequisites: Run after the main scRNA object has been prepared.

rm(list = ls())
library(ggplot2)
library(Seurat)
library(SeuratObject)
library(colorRamp2)
library(ComplexHeatmap)

scRNA <- readRDS("scRNA.rds")
table(scRNA$celltype)
scRNA$celltype <- factor(scRNA$celltype, levels = c("T&NKCell","BCell","PlasmaCell","MyeloidCell",
                                                    "Fibroblast","EndothelialCell","EpithelialCell"))
mean_gene_exp <- AverageExpression(scRNA,
                                   features = c(
                                     "CD2","CD3D",
                                     "GNLY","KLRD1", 
                                     "CD79A","MS4A1",  
                                     "JCHAIN","MZB1",
                                     "LYZ","CD14",  
                                     "COL1A1","LUM",   
                                     "PECAM1","VWF",  
                                     "EPCAM","KRT15"),  
                                   group.by = 'celltype',
                                   slot = 'data') %>%
  data.frame() %>%
  as.matrix()

colnames(mean_gene_exp) <- c("T&NKCell",
                             "BCell","PlasmaCell",
                             "MyeloidCell",
                             "Fibroblast","EndothelialCell",
                             "EpithelialCell")
htdf <- t(scale(t(mean_gene_exp),scale = T,center = T))
col_fun = colorRamp2(c(-2, 0, 2), c("#6DCCFD", "white", "#F484AE"))
column_ha = HeatmapAnnotation(cluster = colnames(htdf),
                              col = list(cluster = col_celltype))
# plot
pdf("Heatmap_celltype_expression.pdf", width = 5, height = 7)
Heatmap(htdf,
        name = "Z-score",
        cluster_columns = F,cluster_rows = F,
        row_names_gp = gpar(fontface = 'italic',fontsize = 10),
        row_names_side = 'left',
        border = T,
        rect_gp = gpar(col = "white", lwd = 1),
        column_names_side = 'top',
        column_names_rot = 45,
        top_annotation = column_ha,
        col = col_fun,
        width = 3, 
        height = 3 
)
dev.off()