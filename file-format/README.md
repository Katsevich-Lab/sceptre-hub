# The sceptre file format: specification and prototype

> **Draft.** It has not been reviewed and is not a contract. The format is
> version 0.1 until sceptre's first release that uses it, which freezes it as
> 1.0.

An on-disk format for Perturb-seq data and analyses: the same files in R and
Python, with one data file and any number of analyses of it.

The format is defined by `SPEC.md` and `spec/`, and the validators decide
whether a file conforms. The code in `format/` is a prototype by one author,
not a reference implementation: its reader and writer are to move into sceptre
and pysceptre, which become the two reference implementations.

| Folder or file | What it is |
|---|---|
| `REVIEW.md` | One page for reviewers: the format in ten lines, what is expensive to change later, the open questions, and what can be skipped |
| `SPEC.md` | The specification: scope, terms, and numbered statements sorted into file rules, behaviour rules and guidance; no reasons |
| `RATIONALE.md` | The reasons: every choice with its alternatives and evidence, comparisons with the IGVF pipeline's format and with pertpy, what is fixed and what can grow, and the measurements |
| `spec/` | The machine-readable parts: `screen.json` (the structure of a data file), `tables.json` (every table's columns, types and key) and JSON Schemas for the JSON records |
| `format/` | The prototype: reader, writer and validators in R (`sceptre_format.R`) and Python (`sceptre_format.py`), with the same names. It uses neither sceptre nor pysceptre |
| `example/` | Files in the format: a high-MOI dataset with six analyses, a low-MOI dataset with one, and a small synthetic data file with every kind of column. `example/*/expected/` holds values that sceptre itself computed for each analysis |
| `conformance/` | Files that conform and files that each break a rule, with what a validator must report (`cases.json`) |
| `verify/` | The verification suite: every check of the format, and the scripts that build and check `conformance/` |
| `notebook.qmd`, `render.sh` | One document that walks through the format in R and in Python, in tabs. Its render stops if a comparison fails |
| `scale-test/` | Recorded results of one run on a real screen of 207,324 cells |

Start with `REVIEW.md`, then `SPEC.md` for what the format is and
`RATIONALE.md` for why.

## Running the checks

```sh
verify/run.sh     # the verification suite; SCEPTRE_FORMAT_PYTHON picks the Python to use
./render.sh       # the notebook
```

Both use standard libraries only and the committed files in `example/` and
`conformance/` as they are.

| Check | What it shows |
|---|---|
| conformance | Both validators and both readers on every file of `conformance/` |
| `encodings` | Every kind of value, written by each language's standard libraries and read by the other's. The kinds the format allows come back correct; the kinds it forbids still break, as recorded in `verify/expected_encodings.json` |
| `two_writers` | R and Python each write the same seeded random content. The files are compared array by array and column by column, and each language reads the other's file and compares it with the source |
| `independent_readers` | Tools outside this code read the files: anndata, scanpy, pandas, polars and duckdb in Python; HDF5Array and arrow in R |
| `wide_index` | The 64-bit branch of the index type (S4), with the limit lowered because no test file reaches it |
| `derived` | What a reader computes from the files (cells per gRNA and per target, cells in use, pairs passing QC, log2 fold change, significance, control groups) against the values in `example/*/expected/` |

They need, in R: `MuData`, `rhdf5`, `HDF5Array`, `SingleCellExperiment`,
`MultiAssayExperiment`, `Matrix`, `jsonlite`, `nanoparquet`, `jsonvalidate`,
`arrow` (`nanoparquet` and `jsonvalidate` may be in `./rlib`), and for the
notebook `reticulate` and quarto. In Python: `mudata`, `anndata`, `h5py`,
`pyarrow`, `pandas`, `scipy`, `jsonschema`, `scanpy`, `polars`, `duckdb`.

`verify/README.md` says how another reader or writer is tested against the
same files.

## What is and is not tested

- **Tested:** the two languages write the same files and read each other's;
  tools outside this code read them; two validators agree on every conformance
  file; values a reader derives match what sceptre computed. One run at scale
  covered the data file.
- **Not tested: whether the files hold everything a real analysis needs.**
  That takes rebuilding a sceptre object from them, and is tested when the
  packages adopt the format. `SPEC.md` marks the rules that have no test.
- **The example files are fixed.** Nothing here regenerates them, so a format
  change converts them in place, with code kept beside the change.
- **One platform.** Everything ran on macOS. The notebook runs R and Python in
  one process and is not known to work on Linux.

Generated and not committed: `rlib/`, `verify/out/`,
`conformance/results-*.json` and `notebook.html`.
