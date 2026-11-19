#!/bin/sh
# The verification suite: every check of the format. It needs neither sceptre
# nor pysceptre. It uses the committed files in example/ and conformance/ as
# they are, standard libraries, and the format code in format/.
#     verify/run.sh            every check
#     verify/run.sh format     the checks that use no writer: the validators and
#                              standard readers on the committed files
#     verify/run.sh writers    the checks of the two writers in format/
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
    "$PY" verify/encodings.py write
    Rscript verify/encodings.R
    "$PY" verify/encodings.py read
    "$PY" verify/independent_readers.py
    Rscript verify/independent_readers.R
    "$PY" verify/derived.py
    Rscript verify/derived.R
fi
if [ "$PART" != format ]; then
    "$PY" verify/random_inputs.py
    Rscript verify/two_writers.R
    "$PY" verify/two_writers.py
    "$PY" verify/wide_index.py
    Rscript verify/wide_index.R
    "$PY" verify/writers.py
    Rscript verify/writers.R
fi
"$PY" verify/summary.py "$PART"
