# Input file specification

This repository does not bundle datasets. Users should supply files with the names and structures expected by each workflow.

## scRNA-seq inputs

- `scRNA.rds`: Seurat object with RNA assay and metadata.
- `CD8TCells.rds`: Seurat subset for CD8 analyses.
- `combined_TCR.rds`: TCR clonotype object aligned to the CD8 subset.
- `cDC.rds`: Seurat subset for dendritic-cell analyses.
- `EpithelialCells_clean.rds`: Seurat subset for NMF.
- `Genes_nmf_w_basis.rds`: NMF basis summary object.
- `out_SCENIC.loom`: SCENIC output file.
- `cellchat_high.rds`, `cellchat_low.rds`: CellChat objects.

## Bulk RNA-seq inputs

- `norm_data_count.csv`: count matrix with genes in rows and samples in columns.
- `group.txt`: tab-delimited metadata with at least `sample`, `group`, and `dataset` columns.
- `Primary_exp.csv`: expression matrix for ssGSEA.
- `signature.csv`: signature gene membership table.
- `survival_result.csv`: survival table with time, status, and stratification variables required by the script.

## Microbiome inputs

- `otu_table_end.csv`: OTU abundance table.
- `tax_table_end.csv`: taxonomy table expected to provide seven taxonomic ranks.
- `metadata.txt`: sample metadata table.
- `physeq.end.rds`: phyloseq object derived from the OTU and metadata tables.
- `Cor_IBL&ThreeGenus_merge.csv`: merged biomarker-correlation input table.

## Spatial inputs

- `cell_anno_insitutype.h5ad`: AnnData object with spatial coordinates and cell annotations.
- `insitutype_anno.csv`: optional annotation table for the Python helper; must contain `cell_barcode` and `cluster_assignment`.
- `cosmx_tx_file.csv`: bacterial coordinate table; must contain `target`, `x_global_px`, and `y_global_px`.
- `high_low_mid_tex_mreg.csv`: summary table for restratification plots.
