#!/bin/sh
# The package-free suite: every check of the format that needs neither sceptre
# nor pysceptre. It uses the committed files in example/ and conformance/ as
# they are, standard libraries, and the format code in format/.
#     verify/run.sh
# Python is taken from SCEPTRE_FORMAT_PYTHON (default: python3).
set -e
cd "$(dirname "$0")/.."
export LANG=en_US.UTF-8
PY="${SCEPTRE_FORMAT_PYTHON:-python3}"
rm -rf verify/out
"$PY" check_conformance.py
Rscript check_conformance.R
"$PY" verify/encodings.py write
Rscript verify/encodings.R
"$PY" verify/encodings.py read
"$PY" verify/random_inputs.py
Rscript verify/two_writers.R
"$PY" verify/two_writers.py
"$PY" verify/independent_readers.py
Rscript verify/independent_readers.R
"$PY" verify/derived.py
Rscript verify/derived.R
"$PY" verify/summary.py
