# Environment manifests

This directory lists the software packages required by each workflow.

## Files

- `r-scrna-packages.txt`
- `r-bulk-rna-packages.txt`
- `r-microbiome-packages.txt`
- `r-spatial-packages.txt`
- `python-spatial-requirements.txt`

## Notes

- The original repository did not include lock files or version-pinned environment manifests.
- Package lists here were derived from the workflow scripts.
- Where version pinning is required for a formal release, maintainers should convert these lists into pinned `renv`, `conda`, or `requirements.txt` snapshots generated from the production environment.
