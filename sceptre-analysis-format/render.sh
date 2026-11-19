#!/bin/sh
# Render overview.qmd to overview.html, the introduction to the format. It
# reports the results of the verification suite, so run verify/run.sh first.
# R and Python run in one process (knitr with reticulate). Tested on macOS only.
# Python is taken from SCEPTRE_FORMAT_PYTHON (default: python3). The data file format is taken from ../perturb-seq-data-format.
#     ./render.sh
set -e
cd "$(dirname "$0")"
export LANG=en_US.UTF-8
export RETICULATE_PYTHON="${SCEPTRE_FORMAT_PYTHON:-$(command -v python3)}"
quarto render overview.qmd --execute-dir .
