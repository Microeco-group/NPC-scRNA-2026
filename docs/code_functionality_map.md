# Code functionality to manuscript mapping

Use this file to map repository workflows and scripts to the manuscript text required by the Nature Research checklist.

| Repository component | Scientific function | Suggested manuscript location |
| --- | --- | --- |
| `workflows/scrna/01_seurat_integration.R` | scRNA-seq preprocessing, integration, clustering | Methods: single-cell RNA-seq analysis |
| `workflows/scrna/02_umap_visualization.R` to `05_celltype_dendrogram.R` | scRNA-seq visualization and annotation summaries | Methods or figure legends |
| `workflows/scrna/06_add_module_score.R` to `15_cellphonedb_visualization.R` | downstream scRNA-seq functional analyses | Methods: downstream single-cell analyses |
| `workflows/bulk-rna/01_combat_seq_batch_correction.R` to `07_ssgsea_scoring.R` | bulk RNA-seq preprocessing and downstream interpretation | Methods: bulk RNA-seq analysis |
| `workflows/microbiome/01_alpha_diversity.R` to `06_correlation_analysis.R` | microbiome diversity, differential taxa, ecological modeling, correlations | Methods: microbiome analysis |
| `workflows/spatial/01_synora_boundary_detection.R` to `04_restratification_comparison.R` | spatial transcriptomics boundary and bacteria-distance analyses | Methods: spatial transcriptomics analysis |
| `env/` and `config/` | environment setup and run configuration | Code availability / software availability statements |
| `docs/reproduction.md` | end-to-end execution guidance | Code availability statement or supplementary methods |

Update the manuscript section labels in the third column before submission so they match the final paper text exactly.
