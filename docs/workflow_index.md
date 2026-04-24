# Workflow index

This document is the master index for workflow purpose, execution order, inputs, outputs, and dependencies.

## 1. scRNA-seq workflow

Directory: `<repository_root>/workflows/scrna`

### Main execution order
1. `01_seurat_integration.R`
2. `02_umap_visualization.R`
3. `03_marker_heatmap.R`
4. `04_roe_analysis.R`
5. `05_celltype_dendrogram.R`
6. `06_add_module_score.R`
7. `07_milor_differential_abundance.R`
8. `08_screpertoire_tcr_analysis.R`
9. `09_monocle_pseudotime.R`
10. `10_cellchat_communication.R`
11. `11_scmetabolism_analysis.R`
12. `12_nmf_program_discovery.R`
13. `13_meta_program_analysis.R`
14. `14_scenic_regulatory_network.R`
15. `15_cellphonedb_visualization.R`

### Entry inputs
- `scRNA.rds`
- `CD8TCells.rds`
- `combined_TCR.rds`
- `cDC.rds`
- `EpithelialCells_clean.rds`
- `Genes_nmf_w_basis.rds`
- `out_SCENIC.loom`
- `cellchat_high.rds`
- `cellchat_low.rds`

### Key dependency notes
- `01_seurat_integration.R` is the primary entry point.
- Steps 02-05 are visualization/annotation-dependent downstream analyses.
- Steps 06-15 consume curated subsets or outputs generated upstream or outside the repository.

See `<repository_root>/docs/scrna_workflow.md`.

## 2. Bulk RNA-seq workflow

Directory: `<repository_root>/workflows/bulk-rna`

### Main execution order
1. `01_combat_seq_batch_correction.R`
2. `02_differential_expression.R`
3. `03_mcp_counter_abundance.R`
4. `04_gsea_pathway_analysis.R`
5. `05_decoupleR_tf_activity.R`
6. `06_survival_analysis.R`
7. `07_ssgsea_scoring.R`

### Entry inputs
- `norm_data_count.csv`
- `group.txt`
- `Primary_exp.csv`
- `signature.csv`
- `survival_result.csv`

### Key dependency notes
- The mandatory chain is 01 -> 02 -> 04/05.
- `03_mcp_counter_abundance.R`, `06_survival_analysis.R`, and `07_ssgsea_scoring.R` are side branches.

See `<repository_root>/docs/bulk_rna_workflow.md`.

## 3. Microbiome workflow

Directory: `<repository_root>/workflows/microbiome`

### Main execution order
1. `01_alpha_diversity.R`
2. `02_beta_diversity.R`
3. `03_lefse_analysis.R`
4. `04_neutral_model_fit.R`
5. `05_cooccurrence_network.R`
6. `06_correlation_analysis.R`

### Entry inputs
- `otu_table_end.csv`
- `tax_table_end.csv`
- `metadata.txt`
- `physeq.end.rds`
- `Cor_IBL&ThreeGenus_merge.csv`

### Key dependency notes
- Alpha/beta diversity establish the core abundance context.
- LEfSe, neutral model fitting, and co-occurrence network analysis can be run as downstream branches.

See `<repository_root>/docs/microbiome_workflow.md`.

## 4. Spatial workflow

Directory: `<repository_root>/workflows/spatial`

### Main execution order
1. `01_synora_boundary_detection.R`
2. `02_cooccurrence_by_boundary_distance.R`
3. `03_fusobacterium_distance.py`
4. `04_restratification_comparison.R`

### Entry inputs
- `cell_anno_insitutype.h5ad`
- `insitutype_anno.csv` (optional for the Python helper)
- `cosmx_tx_file.csv`
- `high_low_mid_tex_mreg.csv`

### Key dependency notes
- `01_synora_boundary_detection.R` must run before `02_cooccurrence_by_boundary_distance.R`.
- `03_fusobacterium_distance.py` is independent, but expects the same spatial coordinate conventions.

See `<repository_root>/docs/spatial_workflow.md`.
