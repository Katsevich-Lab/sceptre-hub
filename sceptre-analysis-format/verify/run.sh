#!/bin/sh
# The verification suite of the sceptre analysis format. It needs neither
# sceptre nor pysceptre. It uses the committed files in example/ and
# conformance/ as they are, standard libraries, the code in format/, and the
# data file format in ../perturb-seq-data-format (its code, examples and
# conformance files).
#     verify/run.sh            every check
#     verify/run.sh format     the checks that use no writer: the validators and
#                              standard readers on the committed files
#     verify/run.sh writers    the checks of the two writers in format/, which
#                              write a data file and an analysis of it
# Python is taken from SCEPTRE_FORMAT_PYTHON (default: python3).
set -e
cd "$(dirname "$0")/.."
export LANG=en_US.UTF-8
PY="${SCEPTRE_FORMAT_PYTHON:-python3}"
PART="${1:-all}"
rm -rf verify/out
if [ "$PART" != writers ]; then
    "$PY" verify/check_conformance.py
    Rscript verify/check_conformance.R
    "$PY" verify/independent_readers.py
    Rscript verify/independent_readers.R
    "$PY" verify/derived.py
    Rscript verify/derived.R
fi
if [ "$PART" != format ]; then
    "$PY" verify/random_inputs.py
    Rscript verify/two_writers.R
    "$PY" verify/two_writers.py
    "$PY" verify/writers.py
    Rscript verify/writers.R
fi
"$PY" verify/summary.py $([ "$PART" = writers ] && echo writers)
