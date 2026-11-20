# The sceptre analysis format

[![Analysis format checks](https://github.com/Katsevich-Lab/sceptre-hub/actions/workflows/sceptre-analysis-format.yml/badge.svg?branch=main)](https://github.com/Katsevich-Lab/sceptre-hub/actions/workflows/sceptre-analysis-format.yml)

> **Draft.** It has not been reviewed. The format is version 0.1 until
> sceptre's first release that uses it, which freezes it as 1.0.

One folder per analysis of a Perturb-seq data file: the settings and tables of
each stage, the same files for R and Python. It builds on the
[Perturb-seq data file format](../perturb-seq-data-format/), which is separate,
has its own version, and knows nothing about analyses.

**Start with [the overview](https://katsevich-lab.github.io/sceptre-hub/sceptre-analysis-format/overview.html).**
It gives the principles and the main choices with their reasons, shows what
an analysis folder holds, how it is read and written, and how the format is checked.
Its source is `overview.qmd`.

Behind it are two reference documents: `SPEC.md` (the rules, numbered) and
`RATIONALE.md` (the full reasons and evidence, under the same numbers).

## Running the checks

```sh
verify/run.sh     # the verification suite; SCEPTRE_FORMAT_PYTHON picks the Python to use
                  # `verify/run.sh format` runs the part that uses no writer, `verify/run.sh writers` the rest
sceptre-suite/run.sh   # the checks that need sceptre 0.99.1: sceptre's analyses, written and rebuilt
```

`verify/run.sh` uses standard libraries only, the committed files in `example/` and
`conformance/` as they are, and the data file format in
`../perturb-seq-data-format` (its code, examples and conformance files). It
needs what that format needs, and in R also `nanoparquet`, `jsonvalidate` and
`arrow` (`nanoparquet` and `jsonvalidate` may be in
`../perturb-seq-data-format/rlib`); in Python also `pyarrow`, `jsonschema`, `polars` and
`duckdb`. `verify/README.md` says how another implementation is tested against
the same files. The suite passes on macOS, Linux and
Windows.
