#!/bin/sh
# Builds the browsable site of the two formats in site/_site, from the
# Markdown files in the two format folders, which are the source: they are
# copied in unchanged and rendered. Needs Quarto. Run from anywhere.
set -e
cd "$(dirname "$0")"
rm -rf _work _site
mkdir _work
cp _quarto.yml index.md rules.lua site.css _work/
for format in perturb-seq-data-format sceptre-analysis-format; do
    mkdir "_work/$format"
    cp "../$format/SPEC.md" "../$format/RATIONALE.md" "../$format/overview.html" "_work/$format/"
done
quarto render _work
mv _work/_site _site
rm -rf _work
