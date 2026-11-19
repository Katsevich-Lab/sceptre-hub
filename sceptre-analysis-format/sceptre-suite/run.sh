#!/bin/sh
# The sceptre suite: the checks of the analysis format that need sceptre.
# It needs sceptre 0.99.1 installed, which is commit 20c33349c of
# Katsevich-Lab/sceptre, and the R libraries of verify/. Run from anywhere.
set -e
cd "$(dirname "$0")/.."
export LANG=en_US.UTF-8
Rscript sceptre-suite/check.R
