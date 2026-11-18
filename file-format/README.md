# The sceptre file format: specification and prototype

> **Draft, not accepted.** This folder is a proposal for a format, with a
> working prototype. It has not been reviewed and is not a contract.
> [`contracts/file-formats.md`](../contracts/file-formats.md), the hub's
> existing format registry, is stale and will be replaced; this draft does not
> follow it.

A proposed on-disk format for Perturb-seq data and analyses, the same files in
R and Python, with one data file and any number of analyses of it.

| Folder or file | What it is |
|---|---|
| `REVIEW.md` | One page for reviewers: the format in ten lines, what is expensive to change later, the open questions, and what can be skipped |
| `SPEC.md` | The specification: scope, terms, and numbered statements sorted into file rules, behaviour rules and guidance; no reasons |
| `RATIONALE.md` | The reasons: every choice with its alternatives, evidence, what it locks in, and how it compares with the IGVF pipeline's format |
| `spec/` | The machine-readable parts: `screen.json` (the structure of a data file), `tables.json` (every table's columns, types and key) and JSON Schemas for the JSON records |
| `format/` | The format code: reader, writer and validators in R (`sceptre_format.R`) and Python (`sceptre_format.py`), with the same names. It uses neither sceptre nor pysceptre |
| `example/` | Files in the format: a high-MOI dataset with seven analyses, a low-MOI dataset with one, and a small synthetic data file with every kind of column. `example/*/expected/` holds what sceptre itself held for each analysis, captured once |
| `conformance/` | Files that conform, files a reader must still accept, and files that each break one rule, with what a validator must report (`cases.json`) |
| `verify/` | **The package-free suite**: every check of the format that needs neither sceptre nor pysceptre |
| `package-suite/` | **The package suite**: checks that build sceptre and pysceptre objects from the files. They depend on unreviewed branch builds of both packages |
| `scale-test/` | The scripts and recorded results of one run on a real screen of 207,324 cells |
| `make_conformance.py`, `check_conformance.R`, `check_conformance.py` | Build `conformance/` from `example/`, and run each language's validators and readers over every case |

Start with `REVIEW.md`. Then `SPEC.md` for what the format is, or sections 1 and 2 of
`RATIONALE.md` for the idea in two pages, then either notebook. In
`RATIONALE.md`:

- section 9 compares the format with IGVF's and with pertpy's conventions;
- section 10 says what is fixed, what can be added later, and which changes to
  the packages would simplify the format;
- section 11 says how the format is formally defined and lists the questions
  that writing it down precisely raised.

The format is a draft (version 0.1) until sceptre's first release that uses it, which freezes it as 1.0.

## Two suites

**The package-free suite** (`verify/run.sh`) is the evidence about the format.
It loads neither sceptre nor pysceptre, uses the committed files in `example/`
and `conformance/` as they are, and runs on a machine that has only standard
libraries:

| Check | What it shows |
|---|---|
| `encodings` | Every kind of value, written by each language's standard libraries and read by the other's. The kinds the format allows come back correct; the kinds it forbids still break, as recorded in `verify/expected_encodings.json` |
| `two_writers` | R and Python each write the same seeded random content. The files are compared array by array and column by column, and each language reads the other's file and compares it with the source |
| `independent_readers` | Tools outside this code read the files: anndata, scanpy, polars and duckdb in Python; HDF5Array and arrow in R |
| conformance | Both validators and both readers on every file of `conformance/` |
| `derived` | What a reader computes from the files (cells per gRNA and per target, cells in use, pairs passing QC, log2 fold change, significance, control groups) against what sceptre held, in `example/*/expected/`. Those values were produced by sceptre once; running the check does not load it |

```sh
verify/run.sh               # SCEPTRE_FORMAT_PYTHON picks the Python to use
```

It needs, in R: `MuData`, `rhdf5`, `HDF5Array`, `SingleCellExperiment`,
`MultiAssayExperiment`, `Matrix`, `jsonlite`, `nanoparquet`, `jsonvalidate`,
`arrow` (`nanoparquet` and `jsonvalidate` may be in `./rlib`). In Python:
`mudata`, `anndata`, `h5py`, `pyarrow`, `pandas`, `scipy`, `jsonschema`,
`scanpy`, `polars`, `duckdb`.

**The package suite** (`package-suite/run.sh`) answers a different question:
do the files hold everything an analysis needs? That can only be tested by
rebuilding each package's object and comparing it with the original. It also
holds checks that are about the packages and not the format: that an analysis
stopped in R and finished in pysceptre gives the same p-values, and that the
two build the same design matrix. Those are parity results.

```sh
package-suite/run.sh                 # rebuild example/ with sceptre, then render both notebooks
package-suite/run.sh --keep-example  # render only
```

It needs the `sceptre3` builds of both packages:

- R sceptre (branch `sceptre3`) from `~/code/research/sceptre-tour-rlib`
  (`SCEPTRE_TOUR_RLIB`);
- pysceptre (fork branch `sceptre3`) from `~/code/research/sceptre-tour-venv`
  (`SCEPTRE_TOUR_PYTHON`);
- quarto, and `reticulate` in R.

The Python notebook renders first, because the R notebook reads what it writes.
This suite is meant to move into the packages when they adopt the format.

Generated and not committed: `rlib/`, `roundtrip/`, `reference-october-layout/`,
`verify/out/`, the rendered notebooks and `package-suite/example-objects.rds`.

## Scope

The format is defined by `SPEC.md` and `spec/`, and the validators decide
whether a file conforms. There is no canonical writer: the standard libraries
do not write conforming files as they come, and any code whose files pass the
validators and read back correctly is a conforming writer. Appendix B of
`RATIONALE.md` lists what was checked and the limits that remain.
