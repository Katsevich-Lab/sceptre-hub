#!/bin/sh
# Builds the browsable site of the two formats in site/_site. The sources are
# the files of the two format folders: SPEC.md and RATIONALE.md, rendered as
# they are, and overview.qmd, which is run. The folders are copied to a work
# folder and rendered there as one Quarto website, so nothing is written into
# them and no HTML of the site is kept in the repository.
#
# The overviews report the results of the verification suites, so run
# verify/run.sh in both folders first. They need what render.sh needs: R, and
# the Python named in SCEPTRE_FORMAT_PYTHON (default: python3). Run from anywhere.
set -e
cd "$(dirname "$0")"
export LANG=en_US.UTF-8
export RETICULATE_PYTHON="${SCEPTRE_FORMAT_PYTHON:-$(command -v python3)}"
rm -rf _work _site
mkdir _work
cp _quarto.yml index.md rules.lua site.css _work/
for format in perturb-seq-data-format sceptre-analysis-format; do
    if [ ! -d "../$format/verify/out/results" ]; then
        echo "Run $format/verify/run.sh first: the overview reports its results." >&2
        exit 1
    fi
    mkdir "_work/$format"
    (cd "../$format" && tar -cf - --exclude ./rlib --exclude ./overview.html --exclude './*_files' --exclude ./.quarto .) | (cd "_work/$format" && tar -xf -)
    [ -d "../$format/rlib" ] && ln -s "$(cd "../$format/rlib" && pwd)" "_work/$format/rlib"
    # a page of a website takes its layout from the site, not from a standalone page
    sed -i.bak -e '/^    toc-location: left$/d' -e '/^    embed-resources: true$/d' "_work/$format/overview.qmd"
    rm "_work/$format/overview.qmd.bak"
done
quarto render _work
mv _work/_site _site
rm -rf _work
