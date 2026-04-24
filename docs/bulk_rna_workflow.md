# Bulk RNA-seq workflow

Directory: `<repository_root>/workflows/bulk-rna`

## Goal

Perform bulk RNA-seq batch correction, differential expression, enrichment analysis, immune abundance estimation, survival analysis, and ssGSEA scoring.

## Required inputs

| File | Purpose |
| --- | --- |
| `norm_data_count.csv` | Raw or normalized count matrix |
| `group.txt` | Sample metadata with sample, group, and dataset columns |
| `Primary_exp.csv` | Expression matrix for ssGSEA |
| `signature.csv` | Signature gene definitions for ssGSEA |
| `survival_result.csv` | Survival analysis input table |

## Execution order

### Main branch
1. `01_combat_seq_batch_correction.R`
2. `02_differential_expression.R`
3. `04_gsea_pathway_analysis.R`
4. `05_decoupleR_tf_activity.R`

### Side branches
- `03_mcp_counter_abundance.R`
- `06_survival_analysis.R`
- `07_ssgsea_scoring.R`

## Main outputs

| Script | Main outputs |
| --- | --- |
| `01_combat_seq_batch_correction.R` | `Bulk_batch_count.csv` |
| `02_differential_expression.R` | `Bulk_batch_differ.csv` |
| `03_mcp_counter_abundance.R` | `MCPcounter_Results.csv`, `box_MCPcounter.pdf` |
| `04_gsea_pathway_analysis.R` | `GSEA_result.pdf` |
| `05_decoupleR_tf_activity.R` | `TF_diff_decoupleR.pdf` |
| `06_survival_analysis.R` | `survival_plot.pdf` |
| `07_ssgsea_scoring.R` | `ssGSEA_result.csv` |

## Common failure points

- `group.txt` missing the `dataset` column required by ComBat-seq.
- Inconsistent sample names between metadata and expression matrices.
- Mismatch between the script's expected differential expression filename (`differ.csv`) and exported output names.
- Missing clinical columns in `survival_result.csv`.

## Figure/output mapping

- Batch-corrected count matrix: `01_combat_seq_batch_correction.R`
- DEG export: `02_differential_expression.R`
- Immune abundance boxplots: `03_mcp_counter_abundance.R`
- Enrichment plots: `04_gsea_pathway_analysis.R`
- TF activity plots: `05_decoupleR_tf_activity.R`
- Survival plot: `06_survival_analysis.R`
