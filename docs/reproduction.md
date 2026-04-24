# Reproduction guide

This document provides a repository-level reproduction order for users who already have the required input files.

## Recommended sequence

1. Prepare environments described in `<repository_root>/env/README.md`.
2. Review file naming and placement rules in `<repository_root>/config/input_file_specification.md`.
3. Run the scRNA-seq core workflow to create or verify the annotated single-cell reference object.
4. Run bulk RNA-seq and microbiome workflows using their corresponding metadata and abundance tables.
5. Run spatial analyses using the spatial `AnnData` object and transcript coordinate tables.
6. Collect generated figures and tables into a manuscript-specific artifact directory under `<repository_root>/results/`.

## Notes

- Several downstream scripts consume subset objects or external tool outputs that are not generated end-to-end inside this repository.
- The repository intentionally does not include demo datasets; reproduction therefore depends on user-supplied data prepared to the documented file specifications.
- If manuscript figure numbering differs from the repository structure, update `<repository_root>/docs/code_functionality_map.md` to reflect the final figure map used for submission.
