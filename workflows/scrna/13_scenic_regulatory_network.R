# Script: 14_scenic_regulatory_network.R
# Purpose: SCENIC regulon analysis and visualization
# Workflow: scRNA
# Required inputs: cDC.rds; out_SCENIC.loom
# Main outputs: TF_specific.pdf; cDC_SCENIC_regulonRSS.csv and related outputs.
# Prerequisites: Requires a completed SCENIC run or an exported loom file.

rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(patchwork)      
library(ComplexHeatmap)
library(circlize)      
library(Seurat)
library(dplyr)
library(SCENIC)
library(AUCell)
library(SCopeLoomR)
library(data.table)


loom <- open_loom('out_SCENIC.loom') 

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)
rownames(regulonAUC)
names(regulons)

seurat.data <- readRDS("cDC.rds")
Idents(seurat.data) <- "celltype3"
DimPlot(seurat.data,reduction = "umap",label=T )

sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]
dim(sub_regulonAUC)
seurat.data
identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(row.names = colnames(seurat.data), 
                           seurat_clusters = as.character(seurat.data$celltype3))
cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        celltype = seurat.data$celltype3)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 

save(sub_regulonAUC,cellTypes,cellClusters,seurat.data,
     file = 'for_rss_and_visual.Rdata')

regulonsToPlot = c("RELA(+)","RELB(+)","NFKB1(+)","NFKB2(+)")
regulonsToPlot %in% row.names(sub_regulonAUC)
seurat.data@meta.data = cbind(seurat.data@meta.data ,t(assay(sub_regulonAUC[regulonsToPlot,])))

mregDC_seurat <- subset(seurat.data, subset = celltype3 == "c41_mregDC_LAMP3")

p1 = DotPlot(seurat.data, features = unique(regulonsToPlot)) + RotatedAxis()
p2 = RidgePlot(seurat.data, features = regulonsToPlot , ncol = 2, cols = c("#eed785","#d4b595","#efe0e7")) 
p3 = VlnPlot(seurat.data, features = regulonsToPlot,pt.size = 0)
p4 = FeaturePlot(seurat.data,features = regulonsToPlot)

wrap_plots(p1,p2,p3,p4)

p2
ggsave("TF_specific.pdf", width = 5, height = 5)

selectedResolution <- "celltype" 
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])

sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
pdf("TFactivity_cDC.pdf", width = 4, height = 15)
Heatmap(
  regulonActivity_byGroup_Scaled, 
  name                         = "z-score",
  col                          = colorRamp2(seq(from=min(regulonActivity_byGroup_Scaled),to=max(regulonActivity_byGroup_Scaled),length=3),rev(c("#e780a6","white","#71c3ee"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
dev.off()


rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 
rss_original <- rss
plotRSS_oneSet(rss, setName = "c41_mregDC_LAMP3")
mat_df <- as.data.frame(rss)

sorted_mat_df <- mat_df[order(-mat_df$c41_mregDC_LAMP3), ]
write.csv(sorted_mat_df,"cDC_SCENIC_regulonRSS.csv")

rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
plotRSS_oneSet(rss, setName = "c41_mregDC_LAMP3")

activity_zscores <- regulonActivity_byGroup_Scaled
head(activity_zscores)
df = do.call(rbind,
             lapply(1:ncol(activity_zscores), function(i){
               dat= data.frame(
                 path  = rownames(activity_zscores),
                 cluster =   colnames(activity_zscores)[i],
                 sd.1 = activity_zscores[,i],
                 sd.2 = apply(activity_zscores[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(10, fc)
rowcn = data.frame(path = top5$cluster) 
n = activity_zscores[top5$path,] 

new_df <- as.matrix(rss_original)

pdf("TFactivity_RSS_original_cDC.pdf", width = 4, height = 15)
Heatmap(
  as.matrix(new_df),  
  name                         = "z-score",
  col                          = colorRamp2(seq(from=min(new_df),to=max(new_df),length=2),rev(c("#e780a6","white"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
dev.off()
