# The Perturb-seq data file format

[![Data file format checks](https://github.com/Katsevich-Lab/sceptre-hub/actions/workflows/perturb-seq-data-format.yml/badge.svg?branch=main)](https://github.com/Katsevich-Lab/sceptre-hub/actions/workflows/perturb-seq-data-format.yml)

> **Draft.** It has not been reviewed. The format is version 0.1 until it is
> frozen as 1.0.

One MuData `.h5mu` file per Perturb-seq screen: counts, cells, features and
gRNAs, the same file for R and Python. It is a set of names and rules on top
of MuData, meant for any software that handles Perturb-seq data. Nothing in it
is specific to one analysis tool.

**Start with [the overview](https://katsevich-lab.github.io/sceptre-hub/file-formats/perturb-seq-data-format/overview.html).**
It gives the principles and the main choices with their reasons, shows what
a data file holds, how it is read and written, and how the format is checked.
Its source is `overview.qmd`.

Behind it are two reference documents: `SPEC.md` (the rules, numbered) and
`RATIONALE.md` (the full reasons and evidence, under the same numbers).

One format that builds on this one is the
[sceptre analysis format](../sceptre-analysis-format/), which stores sceptre's
analyses of a data file.

## Running the checks

```sh
verify/run.sh     # the verification suite; SCEPTRE_FORMAT_PYTHON picks the Python to use
```

It uses standard libraries only and the committed files in `example/` and
`conformance/` as they are. It needs, in R: `MuData`, `rhdf5`, `HDF5Array`, `SingleCellExperiment`, `MultiAssayExperiment`, `Matrix`, `jsonlite` (and for the overview page
`reticulate` and quarto). In Python: `mudata`, `anndata`, `h5py`,
`pandas`, `scipy`, `scanpy`. The suite passes on macOS, Linux and
Windows.
