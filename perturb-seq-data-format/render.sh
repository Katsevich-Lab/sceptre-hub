#!/bin/sh
# Render notebook.qmd to notebook.html, the introduction to the format. It
# reports the results of the verification suite, so run verify/run.sh first.
# R and Python run in one process (knitr with reticulate). Tested on macOS only.
# Python is taken from SCEPTRE_FORMAT_PYTHON (default: python3).
#     ./render.sh
set -e
cd "$(dirname "$0")"
export LANG=en_US.UTF-8
export RETICULATE_PYTHON="${SCEPTRE_FORMAT_PYTHON:-$(command -v python3)}"
quarto render notebook.qmd --execute-dir .
