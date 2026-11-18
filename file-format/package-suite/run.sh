#!/bin/sh
# The package suite: rebuild the examples with sceptre, and render the two
# notebooks, which build sceptre and pysceptre objects from the files. It needs
# the sceptre3 builds of both packages (see README.md). Run from anywhere.
# The Python notebook runs first: the R notebook reads what it writes.
#     package-suite/run.sh                  rebuild example/, then render
#     package-suite/run.sh --keep-example   render only
set -e
cd "$(dirname "$0")/.."
export LANG=en_US.UTF-8
export SCEPTRE_TOUR_RLIB="${SCEPTRE_TOUR_RLIB:-$HOME/code/research/sceptre-tour-rlib}"
export SCEPTRE_TOUR_PYTHON="${SCEPTRE_TOUR_PYTHON:-$HOME/code/research/sceptre-tour-venv/bin/python}"
export RETICULATE_PYTHON="$SCEPTRE_TOUR_PYTHON"
if [ "$1" != "--keep-example" ]; then
    Rscript package-suite/make_example.R
    "$SCEPTRE_TOUR_PYTHON" make_conformance.py
fi
quarto render package-suite/notebook-python.qmd --execute-dir .
quarto render package-suite/notebook-r.qmd --execute-dir .
