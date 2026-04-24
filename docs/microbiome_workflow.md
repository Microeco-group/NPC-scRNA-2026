# Microbiome workflow

Directory: `<repository_root>/workflows/microbiome`

## Goal

Compute diversity metrics, differential taxa, neutral community fits, co-occurrence networks, and biomarker correlations from microbiome abundance data.

## Required inputs

| File | Purpose |
| --- | --- |
| `otu_table_end.csv` | OTU abundance table |
| `tax_table_end.csv` | Taxonomy table |
| `metadata.txt` | Sample metadata |
| `physeq.end.rds` | Prepared phyloseq object |
| `Cor_IBL&ThreeGenus_merge.csv` | Correlation summary table |

## Execution order

### Core branch
1. `01_alpha_diversity.R`
2. `02_beta_diversity.R`
3. `03_lefse_analysis.R`
4. `04_neutral_model_fit.R`
5. `05_cooccurrence_network.R`
6. `06_correlation_analysis.R`

## Main outputs

| Script | Main outputs |
| --- | --- |
| `01_alpha_diversity.R` | `Plot_alpha.pdf` |
| `02_beta_diversity.R` | `Beta_PCoA_uunifrac.pdf` |
| `03_lefse_analysis.R` | `lefse_diff.csv`, `lefse.pdf` |
| `04_neutral_model_fit.R` | `neutral_fit.pdf` |
| `05_cooccurrence_network.R` | `Co_Network_High&Low.pdf` |
| `06_correlation_analysis.R` | `Cor_IBL&ThreeGenus_rel.pdf` |

## Common failure points

- Taxonomy table has fewer than the expected seven columns.
- OTU table and metadata sample names do not match.
- `physeq.end.rds` is missing required components for distance calculations.
- Correlation analysis input table has already been filtered or reshaped differently from the script expectation.

## Figure/output mapping

- Diversity figures: `01_alpha_diversity.R`, `02_beta_diversity.R`
- Differential abundance plots: `03_lefse_analysis.R`
- Neutral model plot: `04_neutral_model_fit.R`
- Network figure: `05_cooccurrence_network.R`
- Correlation plot: `06_correlation_analysis.R`
