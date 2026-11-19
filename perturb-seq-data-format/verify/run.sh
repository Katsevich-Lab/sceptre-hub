#!/bin/sh
# The verification suite of the Perturb-seq data file format. It uses the
# committed files in example/ and conformance/ as they are, standard libraries,
# and the code in format/. It needs no analysis tool.
#     verify/run.sh
# Python is taken from SCEPTRE_FORMAT_PYTHON (default: python3).
set -e
cd "$(dirname "$0")/.."
export LANG=en_US.UTF-8
PY="${SCEPTRE_FORMAT_PYTHON:-python3}"
rm -rf verify/out
"$PY" verify/check_conformance.py
Rscript verify/check_conformance.R
"$PY" verify/encodings.py write
Rscript verify/encodings.R
"$PY" verify/encodings.py read
"$PY" verify/independent_readers.py
Rscript verify/independent_readers.R
"$PY" verify/wide_index.py
Rscript verify/wide_index.R
"$PY" verify/writers.py
Rscript verify/writers.R
"$PY" verify/summary.py
