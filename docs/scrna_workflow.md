# scRNA-seq workflow

Directory: `<repository_root>/workflows/scrna`

## Goal

Process and analyze single-cell RNA-seq data for clustering, visualization, abundance shifts, TCR analysis, cell-cell communication, metabolism, trajectory, NMF, and regulon discovery.

## Required inputs

| File | Purpose |
| --- | --- |
| `scRNA.rds` | Main Seurat object used by the core workflow |
| `CD8TCells.rds` | CD8 T-cell subset for module scoring, miloR, and scRepertoire |
| `combined_TCR.rds` | TCR clonotype input for scRepertoire |
| `cDC.rds` | Dendritic-cell subset for scMetabolism and SCENIC |
| `EpithelialCells_clean.rds` | Epithelial-cell subset for NMF |
| `Genes_nmf_w_basis.rds` | NMF-derived gene basis table |
| `out_SCENIC.loom` | SCENIC result file |
| `cellchat_high.rds` / `cellchat_low.rds` | Group-specific CellChat objects |
| `CD8_SignatureGeneSet.csv` | Signature definitions |
| `row.txt`, `colum.txt` | CellPhoneDB-derived tables |

## Execution order

### Core branch
1. `01_umap_visualization.R`
2. `02_marker_heatmap.R`
3. `03_roe_analysis.R`
4. `04_celltype_dendrogram.R`

### Downstream branch analyses
- `05_add_module_score.R`
- `06_milor_differential_abundance.R`
- `07_screpertoire_tcr_analysis.R`
- `08_monocle_pseudotime.R`
- `09_cellchat_communication.R`
- `10_scmetabolism_analysis.R`
- `11_nmf_program_discovery.R`
- `12_meta_program_analysis.R`
- `13_scenic_regulatory_network.R`
- `14_cellphonedb_visualization.R`

## Main outputs

| Script | Main outputs |
| --- | --- |
| `01_umap_visualization.R` | `Dimplot_scRNA.pdf` |
| `03_roe_analysis.R` | `Roe_result.rds` |
| `06_milor_differential_abundance.R` | `miloR_UMAP.pdf`, `miloR_DA.pdf` |
| `07_screpertoire_tcr_analysis.R` | `cloneType.pdf`, `CD8TCR_occupiedscRepertoire_table.csv` |
| `08_monocle_pseudotime.R` | `Pseudotime_celltype.pdf` |
| `09_cellchat_communication.R` | multiple CellChat comparison plots |
| `10_scmetabolism_analysis.R` | `cDC_scMetabolism_type.rds` |
| `11_nmf_program_discovery.R` | `res_nmf_sampleX.rds` |
| `12_meta_program_analysis.R` | `Robust_NMF_MP.pdf`, `MP_list.rds` |
| `13_scenic_regulatory_network.R` | `TF_specific.pdf`, `cDC_SCENIC_regulonRSS.csv` |
| `14_cellphonedb_visualization.R` | `npc_project_merge.txt`, `npc_project_filt.txt`, `CellphoneDB_bubble.pdf` |

## Common failure points

- Missing subset objects such as `CD8TCells.rds` or `cDC.rds`.
- Missing metadata columns in the Seurat object.
- External tools or outputs not yet prepared (`CellChat`, `SCENIC`, `CellPhoneDB`).
- Online dependency on Ensembl/biomaRt for `11_scmetabolism_analysis.R`.

## Figure/output mapping

- UMAP figure generation: `02_umap_visualization.R`
- Marker heatmaps: `03_marker_heatmap.R`
- Dendrograms: `05_celltype_dendrogram.R`
- TCR figures: `08_screpertoire_tcr_analysis.R`
- Cell-cell communication figures: `10_cellchat_communication.R`
- Regulon figures: `14_scenic_regulatory_network.R`
