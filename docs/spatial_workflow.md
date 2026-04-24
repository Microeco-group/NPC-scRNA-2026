# Spatial workflow

Directory: `<repository_root>/workflows/spatial`

## Goal

Analyze spatial boundary structure, distance-stratified cell abundance, bacteria-to-cell distances, and restratified comparisons.

## Required inputs

| File | Purpose |
| --- | --- |
| `cell_anno_insitutype.h5ad` | Main spatial `AnnData` object |
| `insitutype_anno.csv` | Optional annotation table for the Python helper |
| `cosmx_tx_file.csv` | Bacterial transcript coordinate table |
| `high_low_mid_tex_mreg.csv` | Restratification summary table |

## Execution order

### Main branch
1. `01_synora_boundary_detection.R`
2. `02_cooccurrence_by_boundary_distance.R`

### Independent helper branch
3. `03_fusobacterium_distance.py`
4. `04_restratification_comparison.R`

## Main outputs

| Script | Main outputs |
| --- | --- |
| `01_synora_boundary_detection.R` | `cosmx_synora.csv` |
| `02_cooccurrence_by_boundary_distance.R` | `combined_data_Tex.csv` and distance-stratified plots |
| `03_fusobacterium_distance.py` | `output/tex_distance/disfuso_tex.h5ad`, `output/mreg_distance/disfuso_mreg.h5ad` |
| `04_restratification_comparison.R` | `Proportion_restratification.pdf` |

## Common failure points

- `AnnData` object lacks `CenterX_global_px`, `CenterY_global_px`, or `cluster_assignment` fields.
- Annotation CSV for the Python helper does not contain `cell_barcode` and `cluster_assignment`.
- `cosmx_tx_file.csv` does not include `target`, `x_global_px`, and `y_global_px`.
- No rows match the requested bacterial target or cluster labels.

## Figure/output mapping

- Boundary result table: `01_synora_boundary_detection.R`
- Co-occurrence/distance plots: `02_cooccurrence_by_boundary_distance.R`
- Fusobacterium distance subsets: `03_fusobacterium_distance.py`
- Restratification summary plot: `04_restratification_comparison.R`
