# Script: 09_monocle_pseudotime.R
# Purpose: Monocle trajectory and pseudotime analysis
# Workflow: scRNA
# Required inputs: scRNA.rds with curated variable features and annotations.
# Main outputs: Pseudotime_celltype.pdf and related trajectory plots.
# Prerequisites: Requires a processed Seurat object.

#载入R包
library(Seurat)
library(monocle)
library(dplyr)
library(base)
library(igraph)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(viridis)
library(tidyverse)
library(tidydr)
library(ggforce)
library(ggrastr)

scRNA <- readRDS("scRNA.rds")
class(scRNA)
DimPlot(scRNA)

data <- as(as.matrix(scRNA@assays$RNA@counts),'sparseMatrix')
pd <- new('AnnotatedDataFrame',data = scRNA@meta.data)
FeatureData <- data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new('AnnotatedDataFrame',data = FeatureData)
scRNA_cds <- newCellDataSet(data,
                                phenoData = pd, 
                                featureData = fd,
                                expressionFamily = negbinomial.size())
class(scRNA_cds)

scRNA_cds <- estimateSizeFactors(scRNA_cds)
scRNA_cds <- estimateDispersions(scRNA_cds)
nrow(scRNA_cds)

var.genes <- VariableFeatures(scRNA)
scRNA_cds <- setOrderingFilter(scRNA_cds, var.genes)
plot_ordering_genes(scRNA_cds)

scRNA_cds <- reduceDimension(scRNA_cds, max_components = 2, method = 'DDRTree')
scRNA_cds_temp <- orderCells(scRNA_cds,reverse = F)  

p1 <- plot_cell_trajectory(scRNA_cds_temp, color_by = "Pseudotime", cell_size = 1, show_branch_points = F) +
  scale_colour_gradient(low = "#FAF5A9", high = "#C52F2D") + theme(legend.position = "right") +
  scale_x_reverse()+
  theme(aspect.ratio = 0.5)
p2 <- plot_cell_trajectory(scRNA_cds_temp, color_by = "celltype3",cell_size = 1,show_branch_points=F)+
  scale_color_manual(values= col_celltype3) + theme(legend.position = "right") +
  scale_x_reverse()+
  theme(aspect.ratio = 0.5)
p1/p2 + plot_layout(guides = "collect")
ggsave(file = "Pseudotime_celltype.pdf", width = 6,height = 4)
