# Script: 05_celltype_dendrogram.R
# Purpose: Cell-type dendrogram and circular visualization
# Workflow: scRNA
# Required inputs: Annotated Seurat object with celltype3 labels.
# Main outputs: Dendrogram PDF files written by the script.
# Prerequisites: Requires final cell-type annotations.

rm(list = ls())
library(dplyr)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)
library(dendextend)
library(circlize)
library(ggsci)
library(viridis)

table(scRNA$celltype3)
ave_exp <- AverageExpression(scRNA,
                             group.by = 'celltype3',
                             slot = 'data') %>%
  data.frame() %>%
  as.matrix()
colnames(ave_exp)

result <- sub("^RNA\\.", "", colnames(ave_exp))
print(result)
colnames(ave_exp) <- result
ave_exp <- t(ave_exp)

range(rowSums(ave_exp)) # 10000 10000
ave_exp <- log1p(ave_exp)
ave_exp[1:5, 1:5]
dim(ave_exp) 

# label color
col_celltype3 <- c("#5eb4f3","#ea94cd","#b2df8a","#D9CE99","#f5bcc9","#84d7f7","#815AA8","#80c1c4",   #CD8T
                   "#c5d9df","#e16db7","#908ebc","#af88bb","#7fc28e","#f07590","#dfb9d5","#b099b5","#5394c3","#b8a89f","#917393",  #CD4T
                   "#b589bc","#ecb888","#fc7440","#a032cb","#cd8ab1","#e9d6eb","#76a1ec",   # #Innate_CTL
                   "#f48f89","#f9cabe","#c76da8","#fef6ae","#6276b2","#dedbee","#d9706f","#6f89d3",  #B
                   "#efc1ec","#AB3282","#d8b6f7","#dbebf7",  #Plasma
                   "#eed785","#d4b595","#efe0e7","#eebed6","#f7d39b","#cd7560","#edb869","#e69d8e","#e3edc0","#f8e3db",  #Myeloid
                   "#92CF72","#33a02c","#a3c4a5",   #CAF
                   "#ffff99","#fdbf6f","#f7f4b7",    #EC 
                   "#6a3d9a","#cab2d6")  #Epithelial
colors_to_use <- col_celltype3[1:54]

hc <- as.dist((1 - cor(t(ave_exp))) / 2) %>% # dist(ave_exp) %>%
  hclust(method = "ward.D") %>%
  as.dendrogram() %>%
  set("labels_cex", 0.5)
labels_colors(hc) <- colors_to_use

# Circular dendrogram
pdf("Circular_dendrogram.pdf", width = 4, height = 4)
circlize_dendrogram(hc,
                    labels_track_height = 0.6,
                    dend_track_height = 0.2
)
dev.off()

pdf("Narrow_dendrogram.pdf", width = 4, height = 8) 
par(mar = c(4, 1, 2, 10))
plot(hc, 
     horiz = TRUE,         # 水平方向
     main = "Hierarchical Clustering Dendrogram",
     xlab = "Distance",
     ylab = "",
     axes = FALSE)
dev.off()