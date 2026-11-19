# The sceptre analysis format

[![Analysis format checks](https://github.com/Katsevich-Lab/sceptre-hub/actions/workflows/sceptre-analysis-format.yml/badge.svg?branch=main)](https://github.com/Katsevich-Lab/sceptre-hub/actions/workflows/sceptre-analysis-format.yml)

> **Draft.** It has not been reviewed. The format is version 0.1 until
> sceptre's first release that uses it, which freezes it as 1.0.

One folder per analysis of a Perturb-seq data file: the settings and tables of
each stage, the same files in R and in Python. It builds on the
[Perturb-seq data file format](../perturb-seq-data-format/), which is separate,
has its own version, and knows nothing about analyses.

**Start with [the notebook](https://raw.githack.com/Katsevich-Lab/sceptre-hub/main/sceptre-analysis-format/notebook.html).**
It gives the principles, shows what an analysis folder holds, the choices to
review and their reasons, the questions for reviewers, what has been tested, and what each
folder here holds. Its source is `notebook.qmd`.

Behind it are two reference documents: `SPEC.md` (the rules, numbered) and
`RATIONALE.md` (the full reasons and evidence, under the same numbers).

## Running the checks

```sh
verify/run.sh     # the verification suite; SCEPTRE_FORMAT_PYTHON picks the Python to use
                  # `verify/run.sh format` runs the part that uses no writer, `verify/run.sh writers` the rest
./render.sh       # notebook.html, which reports the suite's results
```

Both use standard libraries only, the committed files in `example/` and
`conformance/` as they are, and the data file format in
`../perturb-seq-data-format` (its code, examples and conformance files). They
need what that format needs, and in R also `nanoparquet`, `jsonvalidate` and
`arrow` (`nanoparquet` and `jsonvalidate` may be in
`../perturb-seq-data-format/rlib`); in Python also `jsonschema`, `polars` and
`duckdb`. `verify/README.md` says how another implementation is tested against
the same files. The suite passes on macOS and on Linux (Ubuntu 24.04).
