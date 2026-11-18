#!/bin/sh
# Render notebook.qmd: one document, R and Python in tabs, standard libraries
# only. R and Python run in one process (knitr with reticulate). Tested on
# macOS only. Python is taken from SCEPTRE_FORMAT_PYTHON (default: python3).
#     ./render.sh
set -e
cd "$(dirname "$0")"
export LANG=en_US.UTF-8
export RETICULATE_PYTHON="${SCEPTRE_FORMAT_PYTHON:-$(command -v python3)}"
quarto render notebook.qmd --execute-dir .
