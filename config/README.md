# Configuration guidance

The workflow scripts in this repository mostly use relative file paths and assume execution from the workflow directory.

## Recommended conventions

1. Change into the relevant workflow directory before running a script.
2. Place required input files directly in that workflow directory unless you update the script paths.
3. Keep generated files in the same workflow directory or in a workflow-local `output/` subdirectory.
4. Use `<repository_root>/config/path-template.env` as a starting point if you want to manage paths with environment variables.
