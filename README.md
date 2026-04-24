# NPC_scRNAseq

This repository contains reproducible analysis scripts for a multi-omics nasopharyngeal carcinoma (NPC) study. The codebase is organized as four workflows:

- scRNA-seq
- bulk RNA-seq
- microbiome
- spatial transcriptomics

The repository has been restructured to align with the Nature Research software checklist **Required content** items, while intentionally excluding bundled demo datasets.

## Repository structure

```text
<repository_root>
├── README.md
├── LICENSE
├── workflows/
│   ├── scrna/
│   ├── bulk-rna/
│   ├── microbiome/
│   └── spatial/
├── docs/
├── env/
├── config/
└── results/
```

See `<repository_root>/docs/workflow_index.md` for the workflow dependency map.

## System requirements

### Operating systems

The scripts are intended for Linux or macOS environments with command-line access to R and Python.

### Language requirements

- Python: the refactored spatial helper script was checked in a Python 3.12 environment.
- R: an R installation is required for the R workflows. The original repository did not record a pinned R version; use a recent Bioconductor-compatible release.

### Package requirements

Per-workflow package lists are provided in `<repository_root>/env/`:

- `<repository_root>/env/r-scrna-packages.txt`
- `<repository_root>/env/r-bulk-rna-packages.txt`
- `<repository_root>/env/r-microbiome-packages.txt`
- `<repository_root>/env/r-spatial-packages.txt`
- `<repository_root>/env/python-spatial-requirements.txt`

### Hardware requirements

No non-standard hardware is required, but the scRNA-seq and spatial workflows may require substantial memory for large `Seurat` or `AnnData` objects.

### Tested environment notes

This sandbox did not include R, so full end-to-end runtime validation of the R workflows was not possible here. The repository now documents the required packages and execution order so maintainers can validate in their target environment.

## Installation guide

### 1. Clone the repository

```bash
git clone https://github.com/LML1qa2ws/NPC_scRNAseq.git
cd NPC_scRNAseq
```

### 2. Prepare software environments

Review `<repository_root>/env/README.md` and install the packages required for the workflow you plan to run.

For Python-dependent spatial analysis:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r <repository_root>/env/python-spatial-requirements.txt
```

For R-dependent workflows, install the packages listed in the corresponding `r-*.txt` file using CRAN/Bioconductor as appropriate.

### 3. Review input and path conventions

Before running any workflow, read:

- `<repository_root>/config/README.md`
- `<repository_root>/config/input_file_specification.md`
- `<repository_root>/config/path-template.env`

### Typical installation time

Environment setup depends on package availability and network speed. A fresh installation for R/Bioconductor-heavy workflows may take substantially longer than a lightweight Python-only setup.

## Demo

This repository does **not** ship a demo dataset. To keep the repository policy-aligned while excluding dataset content, the documented demo is a **dry-run style execution guide** using user-supplied files with the same names and structures described in `<repository_root>/config/input_file_specification.md`.

Minimal demo procedure:

1. Choose one workflow directory under `<repository_root>/workflows/`.
2. Place the required input files in that workflow directory.
3. Run the scripts in the documented order from the workflow directory.
4. Confirm that the documented output files are generated.

Example for the spatial Python helper:

```bash
cd <repository_root>/workflows/spatial
python 03_fusobacterium_distance.py \
  --adata cell_anno_insitutype.h5ad \
  --annotation-csv insitutype_anno.csv \
  --bacteria-csv cosmx_tx_file.csv \
  --output-root output
```

Expected demo outputs:

- `output/tex_distance/disfuso_tex.h5ad`
- `output/mreg_distance/disfuso_mreg.h5ad`

Demo runtime varies with input size and available hardware.

## Instructions for use

### General usage rules

- Run each script from its own workflow directory so relative file paths resolve correctly.
- Keep required input files in the workflow directory unless you modify the scripts to read from alternate paths.
- Use `<repository_root>/docs/workflow_index.md` as the master execution map.

### Workflow entry points

- scRNA-seq: `<repository_root>/workflows/scrna/01_seurat_integration.R`
- bulk RNA-seq: `<repository_root>/workflows/bulk-rna/01_combat_seq_batch_correction.R`
- microbiome: `<repository_root>/workflows/microbiome/01_alpha_diversity.R`
- spatial: `<repository_root>/workflows/spatial/01_synora_boundary_detection.R`

Detailed instructions are provided in:

- `<repository_root>/docs/scrna_workflow.md`
- `<repository_root>/docs/bulk_rna_workflow.md`
- `<repository_root>/docs/microbiome_workflow.md`
- `<repository_root>/docs/spatial_workflow.md`

## Optional reproduction instructions

For a full figure-generation order across workflows, see:

- `<repository_root>/docs/reproduction.md`

## Additional information

### License

This repository now includes `<repository_root>/LICENSE`.

### Source code repository link

- https://github.com/Microeco-group/NPC-scRNA-2026

### Code functionality description location

A code-to-manuscript mapping template is provided at:

- `<repository_root>/docs/code_functionality_map.md`

### Results directory guidance

Expected output locations and recommended artifact organization are described in:

- `<repository_root>/results/README.md`
