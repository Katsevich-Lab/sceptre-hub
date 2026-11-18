#!/bin/sh
# Rebuild the example and the conformance files, check both validators against
# them, and render both notebooks. Run from this folder.
# The Python notebook runs first: the R notebook reads what it writes.
set -e
export LANG=en_US.UTF-8
export SCEPTRE_TOUR_RLIB="${SCEPTRE_TOUR_RLIB:-$HOME/code/research/sceptre-tour-rlib}"
export SCEPTRE_TOUR_PYTHON="${SCEPTRE_TOUR_PYTHON:-$HOME/code/research/sceptre-tour-venv/bin/python}"
export RETICULATE_PYTHON="$SCEPTRE_TOUR_PYTHON"
if [ "$1" != "--keep-example" ]; then
    Rscript make_example.R
fi
"$SCEPTRE_TOUR_PYTHON" make_conformance.py
"$SCEPTRE_TOUR_PYTHON" check_conformance.py
Rscript check_conformance.R
quarto render notebook-python.qmd
quarto render notebook-r.qmd
