# Script: 07_milor_differential_abundance.R
# Purpose: miloR differential abundance analysis
# Workflow: scRNA
# Required inputs: CD8TCells.rds with sample and group metadata.
# Main outputs: miloR_UMAP.pdf; miloR_DA.pdf
# Prerequisites: Requires a prepared CD8 T-cell subset.

rm(list = ls())
library(miloR)
library(SingleCellExperiment)
library(scater)
library(patchwork)
library(dplyr)

CD8TCells <- readRDS("CD8TCells.rds")
CD8T_sce <- as.SingleCellExperiment(CD8TCells)

CD8T_milo <- Milo(CD8T_sce)
CD8T_milo <- buildGraph(CD8T_milo, k = 25, d = 25, reduced.dim = "PCA")
CD8T_milo <- makeNhoods(CD8T_milo, prop = 0.25, k = 25, d=25,
                        refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(CD8T_milo)

CD8T_milo <- countCells(CD8T_milo, meta.data = as.data.frame(colData(CD8T_milo)), sample="sample")
head(nhoodCounts(CD8T_milo))

CD8T_design <- data.frame(colData(CD8T_milo))[,c("sample", "group")]
CD8T_design$group <- factor(CD8T_design$group, levels = c("Low","High")) 
CD8T_design <- distinct(CD8T_design)
rownames(CD8T_design) <- CD8T_design$sample
CD8T_design

CD8T_milo <- calcNhoodDistance(CD8T_milo, d=25, reduced.dim = "PCA")

da_results <- testNhoods(CD8T_milo, design = ~group, design.df = CD8T_design)
head(da_results)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)

CD8T_milo <- buildNhoodGraph(CD8T_milo)
umap_pl <- plotReducedDim(CD8T_milo, dimred = "UMAP", colour_by="group", text_by = "celltype3", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")
umap_pl

nh_graph_pl <- plotNhoodGraphDA(CD8T_milo, da_results, layout="UMAP",alpha=0.1) +
  scale_fill_gradient2(low="#6dccfd",               
                       mid="white",                     
                       high="#fd9aa0",                      
                       name="log2FC")+
  scale_x_continuous(breaks = seq(floor(min(mydata$UMAP_1)), ceiling(max(mydata$UMAP_1)), length.out = 4),expand = c(0.03, 0.03)) + 
  scale_y_continuous(breaks = seq(floor(min(mydata$UMAP_2)), ceiling(max(mydata$UMAP_2)), length.out = 4),expand = c(0.03, 0.03)) + 
  scale_size_continuous(range = c(3,7))+
  coord_fixed(ratio = 1) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.key.height = unit(1,'cm'),
        legend.key.width = unit(0.5,'cm'),
        aspect.ratio = 1) + 
  xlab("")+ylab("")+
  scale_edge_width_continuous(range = c(0.1, 0.5))
nh_graph_pl
ggsave(nh_graph_pl, file="miloR_UMAP.pdf", width=10, height=10)

da_results <- annotateNhoods(CD8T_milo, da_results, coldata_col = "celltype3")
head(da_results)
da_results$celltype3 <- factor(da_results$celltype3,levels = rev(c("c01_CD8_Tnaive_LEF1",
                                                                   "c02_CD8_Tcm_GPR183",
                                                                   "c03_CD8_Trm_CXCR6",
                                                                   "c04_CD8_Tstr_HSPA1A",
                                                                   "c05_CD8_Tisg_ISG15",
                                                                   "c06_CD8_Teff_GZMA",
                                                                   "c07_CD8_Tex_PDCD1",
                                                                   "c08_CD8_Tcycling_MKI67")))

plotDAbeeswarm(da_results, group.by = "celltype3",alpha = 0.1)+
  scale_color_gradient2(midpoint=0, low="#6dccfd", mid="white",high="#fd9aa0", space ="Lab" )
ggsave(file="miloR_DA.pdf", width=7, height=7)
