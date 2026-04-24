# Script: 08_screpertoire_tcr_analysis.R
# Purpose: TCR clonotype analysis with scRepertoire
# Workflow: scRNA
# Required inputs: combined_TCR.rds; CD8TCells.rds
# Main outputs: cloneType.pdf; CD8TCR_occupiedscRepertoire_table.csv; clonotype diversity plots.
# Prerequisites: Requires merged TCR calls and a matching Seurat subset.

rm(list = ls())
library(scRepertoire)
library(Seurat)
library(ggsci)
library(ggplot2)

combined <- readRDS("combined_TCR.rds")
CD8TCells <- readRDS("CD8TCells.rds")
seurat <- combineExpression(combined, CD8TCells,       
                            cloneCall="strict", 
                            proportion = FALSE,        
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
head(seurat)

meta <- as.data.frame(seurat@meta.data)
meta$cloneType[is.na(meta$cloneType)] <- "None"
seurat@meta.data <- meta
seurat$cloneType <- factor(seurat$cloneType, levels = c("None","Single (0 < X <= 1)","Small (1 < X <= 5)",
                                                        "Medium (5 < X <= 20)","Large (20 < X <= 100)","Hyperexpanded (100 < X <= 500)"))
DimPlot(seurat, group.by = "cloneType",label = F,
                cols = c("#eef8fc","#c4d2e5","#a5bbd8","#8e97c2","#7f58a3","#522966"))+
  theme(plot.title = element_blank(),
        aspect.ratio = 1)
ggsave("cloneType.pdf", width = 8, height = 6)

clonoType_cols<- c("#c4d2e5","#a5bbd8","#8e97c2","#7f58a3","#522966")
occupiedscRepertoire(seurat, 
                     x.axis = "celltype3",
                     na.include = T,
                     proportion = T,
                     facet.by = "group")+
  scale_fill_manual(values= rev(clonoType_cols))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table <- occupiedscRepertoire(seurat, 
                              x.axis = "celltype3",
                              proportion = T,
                              facet.by = "group",
                              exportTable = T)
write.csv(table,"CD8TCR_occupiedscRepertoire_table.csv")


table <- clonalDiversity(seurat,       
                            cloneCall = "aa",   
                            n.boots = 100,
                            exportTable = TRUE)
library(ggpubr)
if(T){
  mytheme <- theme(
    axis.title = element_text(size = 12,color ="black"), 
    axis.text = element_text(size= 12,color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1 ),
    legend.text = element_text(size= 12),
    legend.title= element_text(size= 12)
  ) }
ggplot(table, aes(x = group, y = Shannon)) + 
  geom_boxplot(aes(fill = group, color = group),  
               position = position_dodge(0.5), 
               width = 0.5, 
               outlier.alpha = 0) + 
  scale_fill_manual(values = scales::alpha(c("#F484AE", "#6DCCFD"), alpha = 0.6)) + 
  scale_color_manual(values = c("#F484AE", "#6DCCFD")) +  
  theme_bw() + 
  mytheme + 
  stat_compare_means(aes(group = group),
                     label = "p.format",
                     method = "wilcox.test",
                     hide.ns = FALSE) + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm")) + 
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.4,  
    dodge.width = 0.85   
  ), 
  aes(fill = group), 
  color = "white",   
  size = 2,         
  shape = 21,       
  stroke = 0)+
  scale_y_continuous(expand=c(0.1,0.1))
ggsave("clonoType_diversity.pdf",width = 5, height = 4)