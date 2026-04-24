# Script: 02_umap_visualization.R
# Purpose: UMAP visualization for annotated scRNA-seq cells
# Workflow: scRNA
# Required inputs: scRNA.rds with UMAP reduction and metadata columns sample/group/celltype.
# Main outputs: Dimplot_scRNA.pdf
# Prerequisites: Requires a processed Seurat object from the integration workflow.

rm(list = ls())
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(SeuratObject)
library(ggpubr)

scRNA <- readRDS("scRNA.rds")
data1 <- scRNA@meta.data 

data2 <- scRNA@reductions[["umap"]]@cell.embeddings%>%as.data.frame()
mydata <- merge(data2,data1[,c("sample","group","celltype")],by=0,all.x = T) %>% column_to_rownames("Row.names")

col_celltype <- c("#1F78B4","#e0c16f","#bdbcdc","#70b8a7",
                   "#FB9A99","#8ebfd1","#ee9ab9")
names(col_celltype) <- c("BCell","EndothelialCell","EpithelialCell","Fibroblast",
                          "MyeloidCell","PlasmaCell","T&NKCell")

p <- ggplot(data = mydata,aes(UMAP_1,UMAP_2,fill=celltype,colour=celltype)) +
  geom_point(shape=21,size=0.01,alpha=1)+
  scale_fill_manual(values = col_celltype)+  
  scale_colour_manual(values = col_celltype)+
  scale_x_continuous(breaks = seq(floor(min(mydata$UMAP_1)), ceiling(max(mydata$UMAP_1)), length.out = 4),expand = c(0.03, 0.03)) + 
  scale_y_continuous(breaks = seq(floor(min(mydata$UMAP_2)), ceiling(max(mydata$UMAP_2)), length.out = 4),expand = c(0.03, 0.03)) + 
  theme_bw(base_rect_size=1 , base_line_size = 0)+  
  guides(fill=guide_legend(override.aes = list(size=3,alpha=1), ncol =1))+  
  coord_fixed(ratio = 1) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.key.height = unit(1,'cm'),
        legend.key.width = unit(0.5,'cm'))
ggsave(p, file="Dimplot_scRNA.pdf", width=8, height=8)
