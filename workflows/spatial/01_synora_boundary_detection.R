# Script: 01_synora_boundary_detection.R
# Purpose: Boundary detection with Synora
# Workflow: spatial
# Required inputs: cell_anno_insitutype.h5ad
# Main outputs: cosmx_synora.csv
# Prerequisites: Requires an annotated AnnData object with spatial coordinates and cell-type labels.

rm(list = ls())

library(Synora)
library(reticulate)
library(dplyr)

ad <- import("anndata")
adata <- ad$read_h5ad("cell_anno_insitutype.h5ad")
meta <- adata$obs
df <- data.frame(
  Cell_ID = meta$cell_id,
  X = meta$CenterX_global_px,
  Y = meta$CenterY_global_px,
  CT = ifelse(meta$celltype1 == "EpithelialCell", 1, 0)
)

res_boundary <- Synora::GetBoundary(
  INPUT = df,
  CELL_ID_COLUMN = "Cell_ID",
  X_POSITION = "X",
  Y_POSITION = "Y",
  ANNO_COLUMN = "CT",
  RADIUS = 166,
  NEST_SPECIFICITY = 0.25,
  BOUNDARY_SPECIFICITY = 0.05
)

res_dist <- Synora::GetDist2Boundary(
  INPUT = res_boundary,
  CELL_ID_COLUMN = "Cell_ID",
  X_POSITION = "X",
  Y_POSITION = "Y",
  ANNO_COLUMN = "SynoraAnnotation",
  ANNO_OF_BOUNDARY = "Boundary"
)

write.csv(res_dist,"cosmx_synora.csv")

