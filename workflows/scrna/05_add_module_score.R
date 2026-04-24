# Script: 06_add_module_score.R
# Purpose: Signature scoring for CD8 T-cell states
# Workflow: scRNA
# Required inputs: CD8TCells.rds; CD8_SignatureGeneSet.csv
# Main outputs: Signature score heatmap files written by the script.
# Prerequisites: Requires a CD8 T-cell subset extracted from the main Seurat object.

rm(list = ls())
library(readr)
library(scales)
library(pheatmap)
library(dplyr)
library(tidydr)
library(Seurat)
CD8_SignatureGeneSet <- read_csv("CD8_SignatureGeneSet.csv")
DefaultAssay(CD8TCells) <- "RNA"
CD8T_score <- CD8TCells
CD8T_score$celltype3 <- as.factor(CD8T_score$celltype3)
Idents(CD8T_score) <- CD8T_score$celltype3
for (i in 1:20) {
  FunctionName <- colnames(CD8_SignatureGeneSet)[i]
  Feature <- list(CD8_SignatureGeneSet[,i])
  Feature <- lapply(Feature, function(x) x[!is.na(x)])
  CD8T_score <-  AddModuleScore(CD8T_score,
                                features = Feature,
                                ctrl = 5,
                                name = FunctionName)
}
colnames(CD8T_score@meta.data)
colnames(CD8T_score@meta.data)[30:49] <- c("Naive", "Activation_Effector_function", "Exhaustion","Proliferation",
                                           "TCR_Signaling", "Cytotoxicity", "Cytokine_Cytokine_receptor",
                                           "Chemokine_Chemokine receptor", "Senescence", "Anergy",
                                           "NFkB_Signaling", "Stress_response", "MAPK_Signaling", "Adhesion",
                                           "IFN_Response","Oxidative_phosphorylation", "Glycolysis", "Fatty_acid_metabolism",
                                           "Pro_apoptosis", "Anti_apoptosis")
Differentiation <- c("Naive", "Activation_Effector_function", "Exhaustion","Proliferation")
Function <- c("TCR_Signaling", "Cytotoxicity", "Cytokine_Cytokine_receptor",
              "Chemokine_Chemokine receptor", "Senescence", "Anergy",
              "NFkB_Signaling", "Stress_response", "MAPK_Signaling", "Adhesion",
              "IFN_Response")
Metabolism <- c("Oxidative_phosphorylation", "Glycolysis", "Fatty_acid_metabolism")
Apoptosis <- c("Pro_apoptosis", "Anti_apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(CD8T_score$celltype3)),
                              nrow = 20)
colnames(FunctionScoreMatrix) <- levels(CD8T_score$celltype3)
rownames(FunctionScoreMatrix) <- MarkerNameVector

for(ci in 1:ncol(FunctionScoreMatrix)){
  for(ri in 1:nrow(FunctionScoreMatrix)){
    FunctionVec <- as_tibble(CD8T_score@meta.data) %>% pull(MarkerNameVector[ri])
    fv <- mean(FunctionVec[CD8T_score$celltype3 == levels(CD8T_score$celltype3)[ci]])
    FunctionScoreMatrix[ri, ci] <- fv
  }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(
  colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
signatureType_row <- data.frame(Signature.type = c(
  rep("Differentiation", length(Differentiation)),
  rep("Function", length(Function)),
  rep("Metabolism", length(Metabolism)),
  rep("Apoptosis", length(Apoptosis))))
rownames(signatureType_row) <- MarkerNameVector
pdf("Heatmap_Signature_CD8T_celltype3.pdf",width = 7, height = 7)
pheatmap(FunctionScoreMatrix,
         show_colnames = T,
         show_rownames = T,
         annotation_row = signatureType_row,
         gaps_row = c(4, 15, 18),
         cluster_rows = F,
         cluster_cols = F,
         breaks = my.breaks,
         color = my.colors,
         border_color = "NA",
         fontsize = 8,
         cellheight = 15, 
         cellwidth = 15,
         width = 4,
         height = 4)
dev.off()