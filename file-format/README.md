# Two file formats for Perturb-seq

[![File format checks](https://github.com/Katsevich-Lab/sceptre-hub/actions/workflows/file-format.yml/badge.svg?branch=main)](https://github.com/Katsevich-Lab/sceptre-hub/actions/workflows/file-format.yml)

> **Draft.** Neither has been reviewed. Each is at version 0.1 until it is
> frozen as 1.0.

- **The Perturb-seq data file format.** One MuData `.h5mu` file per screen:
  counts, cells, features and gRNAs. It is meant for any software that handles
  Perturb-seq data, and nothing in it is specific to sceptre.
- **The sceptre analysis format.** One folder per analysis of a data file:
  the settings and tables of each stage.

Both are the same files in R and in Python. They have separate versions, and
the data file knows nothing about analyses.

**Start with [the notebook](https://raw.githack.com/Katsevich-Lab/sceptre-hub/main/file-format/notebook.html).**
It gives the principles, shows what is on disk, the choices to review and
their reasons, the questions for reviewers, what has been tested, and what
each folder here holds. Its source is `notebook.qmd`.

Behind it are two reference documents: `SPEC.md` (the rules, numbered) and
`RATIONALE.md` (the full reasons and evidence, under the same numbers).

## Running the checks

```sh
verify/run.sh     # the verification suite; SCEPTRE_FORMAT_PYTHON picks the Python to use
                  # `verify/run.sh format` runs the part that uses no writer, `verify/run.sh writers` the rest
./render.sh       # notebook.html, which reports the suite's results
```

Both use standard libraries only and the committed files in `example/` and
`conformance/` as they are. They need, in R: `MuData`, `rhdf5`, `HDF5Array`,
`SingleCellExperiment`, `MultiAssayExperiment`, `Matrix`, `jsonlite`,
`nanoparquet`, `jsonvalidate`, `arrow` (`nanoparquet` and `jsonvalidate` may
be in `./rlib`), and for the notebook `reticulate` and quarto. In Python:
`mudata`, `anndata`, `h5py`, `pyarrow`, `pandas`, `scipy`, `jsonschema`,
`scanpy`, `polars`, `duckdb`. The suite passes on macOS and on Linux (Ubuntu 24.04).

`verify/README.md` says how another implementation is tested against the same
files.
