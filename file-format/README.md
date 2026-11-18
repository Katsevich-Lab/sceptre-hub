# The sceptre file format: specification and prototype

> **Draft, not accepted.** This folder is a proposal for a format, with a
> working prototype. It has not been reviewed and is not a contract.
> [`contracts/file-formats.md`](../contracts/file-formats.md), the hub's
> existing format registry, is stale and will be replaced; this draft does not
> follow it.

A proposed on-disk format for Perturb-seq data and analyses, the same files in
R and Python, with one data file and any number of analyses of it.

| File | What it is |
|---|---|
| `REVIEW.md` | One page for reviewers: the format in ten lines, what is expensive to change later, the open questions, and what can be skipped |
| `SPEC.md` | The specification: scope, terms, and numbered statements sorted into file rules, behaviour rules and guidance; no reasons |
| `RATIONALE.md` | The reasons: every choice with its alternatives, evidence, what it locks in, and how it compares with the IGVF pipeline's format |
| `spec/` | The machine-readable parts: `screen.json` (the structure of a data file), `tables.json` (every table's columns, types and key) and JSON Schemas for `analysis.json` and the stage records |
| `conformance/` | Files that conform, files a reader must still accept, and files that each break one rule, with what a validator must report (`cases.json`) |
| `notebook-r.html` | The R notebook, rendered: reads the examples and checks every part of the format (`notebook-r.qmd` is its source) |
| `notebook-python.html` | The Python notebook, rendered: the same, section by section (`notebook-python.qmd`) |
| `example/` | Files in the format: a high-MOI dataset with seven analyses, a low-MOI dataset with one, and a small synthetic data file with every kind of column |
| `format/sceptre_format.R` | R reader, writer and validators (`validate_screen()`, `validate_analysis()`) |
| `format/sceptre_format.py` | Python reader, writer and validators, with the same names |
| `make_example.R` | Runs the analyses in R and writes `example/` |
| `make_conformance.py` | Builds `conformance/` from `example/` |
| `check_conformance.R`, `check_conformance.py` | Run each language's validators over every conformance case |
| `scale-test/` | The scripts and recorded results of one run on a real screen of 207,324 cells |
| `render.sh` | Rebuilds the examples and conformance files, checks both validators, and renders both notebooks |

Start with `REVIEW.md`. Then `SPEC.md` for what the format is, or sections 1 and 2 of
`RATIONALE.md` for the idea in two pages, then either notebook. In
`RATIONALE.md`:

- section 9 compares the format with IGVF's and with pertpy's conventions;
- section 10 says what is fixed, what can be added later, and which changes to
  the packages would simplify the format;
- section 11 says how the format is formally defined and lists the questions
  that writing it down precisely raised.

The format is a draft (version 0.1) until sceptre's first release that uses it, which freezes it as 1.0.

Generated and not part of the proposal: `example-objects.rds` (the R objects
the examples were written from), `reference-october-layout/` (the same analyses
in the October layout, an independent reference for the Python checks),
`roundtrip/` (copies and variants written by the notebooks) and `rlib/`
(`nanoparquet`).

## Rebuild

```sh
./render.sh                 # rebuild example/ and conformance/, check, render both notebooks
./render.sh --keep-example  # keep example/; rebuild conformance/, check, render
```

It uses the environments of the `sceptre3` data-structures tour:

- R sceptre (branch `sceptre3`) from `~/code/research/sceptre-tour-rlib`
  (`SCEPTRE_TOUR_RLIB`);
- pysceptre (fork branch `sceptre3`) with `mudata`, `anndata` and `pyarrow`
  from `~/code/research/sceptre-tour-venv` (`SCEPTRE_TOUR_PYTHON`);
- the default R library for `MuData`, `rhdf5`, `jsonlite`, `reticulate`;
- `nanoparquet` and `jsonvalidate` in `./rlib`
  (`install.packages(c("nanoparquet", "jsonvalidate"), lib = "rlib")`);
- `jsonschema` in the Python environment.

The Python notebook renders first, because the R notebook reads what it writes
(`roundtrip/python-copy`, `roundtrip/variants-python` and the analysis
`example/highmoi/analyses/continued_in_python`).

## Scope

The reader/writer code lives here, outside both packages, and builds each
package's object through its internal accessors, on unreviewed branch builds.
The packages are unchanged: where the file's names differ from theirs (the
`rna` modality and its covariates, the gRNA table's `target_` coordinates), the
code here translates.
Appendix B of `RATIONALE.md` lists what was checked and the limits that remain.
