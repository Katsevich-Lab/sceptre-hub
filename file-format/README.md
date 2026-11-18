# The sceptre file format: specification and prototype

> **Draft, not accepted.** This folder is a proposal for a format, with a
> working prototype. It has not been reviewed and is not a contract.

A proposed on-disk format for Perturb-seq data and analyses, the same files in
R and Python, with one data file and any number of analyses of it.

**The code here is a prototype by one author, not a reference implementation.**
It exists to show that the format works in both languages while the format is
still being decided. When that is done, the reader and writer move into
sceptre and pysceptre, which become the two reference implementations, and the
copy here is retired. What stays is the specification, `spec/`, the
conformance files, the validators and a harness for testing any implementation
(`verify/README.md`).

| Folder or file | What it is |
|---|---|
| `REVIEW.md` | One page for reviewers: the format in ten lines, what is expensive to change later, the open questions, and what can be skipped |
| `SPEC.md` | The specification: scope, terms, and numbered statements sorted into file rules, behaviour rules and guidance; no reasons |
| `RATIONALE.md` | The reasons: every choice with its alternatives, evidence, what it locks in, and how it compares with the IGVF pipeline's format |
| `spec/` | The machine-readable parts: `screen.json` (the structure of a data file), `tables.json` (every table's columns, types and key) and JSON Schemas for the JSON records |
| `format/` | The prototype's format code: reader, writer and validators in R (`sceptre_format.R`) and Python (`sceptre_format.py`), with the same names. It uses neither sceptre nor pysceptre. Temporary: the reader and writer are to move into the packages |
| `example/` | Files in the format: a high-MOI dataset with seven analyses, a low-MOI dataset with one, and a small synthetic data file with every kind of column. `example/*/expected/` holds what sceptre itself held for each analysis, captured once. The example files are frozen: the code that made them needed sceptre and was removed after commit `7e10af4` |
| `conformance/` | Files that conform and files that each break a rule, with what a validator must report (`cases.json`). Five of the broken ones store counts or index arrays in another number type, which a reader must still read |
| `verify/` | **The package-free suite**: every check of the format that needs neither sceptre nor pysceptre |
| `scale-test/` | The scripts and recorded results of one run on a real screen of 207,324 cells |
| `make_conformance.py`, `check_conformance.R`, `check_conformance.py` | Build `conformance/` from `example/`, and run each language's validators and readers over every case |

Start with `REVIEW.md`. Then `SPEC.md` for what the format is, or sections 1 and 2 of
`RATIONALE.md` for the idea in two pages. In
`RATIONALE.md`:

- section 9 compares the format with IGVF's and with pertpy's conventions;
- section 10 says what is fixed, what can be added later, and which changes to
  the packages would simplify the format;
- section 11 says how the format is formally defined and lists the questions
  that writing it down precisely raised.

The format is a draft (version 0.1) until sceptre's first release that uses it, which freezes it as 1.0.

## The package-free suite, and what was removed

**The package-free suite** (`verify/run.sh`) is the evidence about the format.
It loads neither sceptre nor pysceptre, uses the committed files in `example/`
and `conformance/` as they are, and runs on a machine that has only standard
libraries:

| Check | What it shows |
|---|---|
| `encodings` | Every kind of value, written by each language's standard libraries and read by the other's. The kinds the format allows come back correct; the kinds it forbids still break, as recorded in `verify/expected_encodings.json` |
| `two_writers` | R and Python each write the same seeded random content. The files are compared array by array and column by column, and each language reads the other's file and compares it with the source |
| `independent_readers` | Tools outside this code read the files: anndata, scanpy, polars and duckdb in Python; HDF5Array and arrow in R |
| `wide_index` | The 64-bit branch of the index type (S4), with the limit lowered because no test file reaches it: both writers, both validators, and R's readers (in memory, MuData's backed reader, HDF5Array) on a file whose index arrays are both 64-bit |
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

**The package suite was removed for now** (Eugene, 18 November). It rebuilt
sceptre and pysceptre objects from the files, and it needed unreviewed
`sceptre3` builds of both packages that exist on one machine, with code that
reached into internals the rebuild of the packages will change. It is in git
history up to commit `7e10af4` and comes back inside the packages when they adopt
the format. What that costs:

- **The sufficiency check is gone.** Rebuilding a real sceptre object is the
  only test that the files hold everything an analysis needs, and it found most
  of the gaps in earlier drafts. The format passed it at every step up to
  commit `7e10af4`. Changes after that commit are not covered until the packages
  test them.
- **The example files are frozen.** The code that generated them ran sceptre.
  A format change from now on converts the committed files in place with
  package-free code, written as part of that change and kept: it is a direct
  test that old files can be carried forward.
- **The frozen expectations stay.** `example/*/expected/` holds what sceptre
  itself held for each analysis, produced once by the removed code.
- **The parity checks are gone**: an analysis stopped in R and finished in
  pysceptre, and the agreement of the two design matrices. Those compare the
  packages, which is not the format's job.

Generated and not committed: `rlib/`, `verify/out/` and
`conformance/results-*.json`.

## Scope

The format is defined by `SPEC.md` and `spec/`, and the validators decide
whether a file conforms. There is no canonical writer: the standard libraries
do not write conforming files as they come, and any code whose files pass the
validators and read back correctly is a conforming writer. Appendix B of
`RATIONALE.md` lists what was checked and the limits that remain.
