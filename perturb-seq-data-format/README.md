# The Perturb-seq data file format

[![Data file format checks](https://github.com/Katsevich-Lab/sceptre-hub/actions/workflows/perturb-seq-data-format.yml/badge.svg?branch=main)](https://github.com/Katsevich-Lab/sceptre-hub/actions/workflows/perturb-seq-data-format.yml)

> **Draft.** It has not been reviewed. The format is version 0.1 until it is
> frozen as 1.0.

One MuData `.h5mu` file per Perturb-seq screen: counts, cells, features and
gRNAs, the same file in R and in Python. It is a set of names and rules on top
of MuData, meant for any software that handles Perturb-seq data. Nothing in it
is specific to one analysis tool.

**Start with [the notebook](https://raw.githack.com/Katsevich-Lab/sceptre-hub/main/perturb-seq-data-format/notebook.html).**
It gives the principles, shows what a data file holds, the choices to review
and their reasons, the questions for reviewers, what has been tested, and
what each folder here holds. Its source is `notebook.qmd`.

Behind it are two reference documents: `SPEC.md` (the rules, numbered) and
`RATIONALE.md` (the full reasons and evidence, under the same numbers).

One format that builds on this one is the
[sceptre analysis format](../sceptre-analysis-format/), which stores sceptre's
analyses of a data file.

## Running the checks

```sh
verify/run.sh     # the verification suite; SCEPTRE_FORMAT_PYTHON picks the Python to use
./render.sh       # notebook.html, which reports the suite's results
```

Both use standard libraries only and the committed files in `example/` and
`conformance/` as they are. They need, in R: `MuData`, `rhdf5`, `HDF5Array`, `SingleCellExperiment`, `MultiAssayExperiment`, `Matrix`, `jsonlite` (and for the notebook
`reticulate` and quarto). In Python: `mudata`, `anndata`, `h5py`, `pyarrow`,
`pandas`, `scipy`, `scanpy`. The suite passes on macOS, Linux and
Windows.
