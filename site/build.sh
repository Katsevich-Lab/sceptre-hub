#!/bin/sh
# Builds the SCEPTRE site in site/_site: a home page, the roadmap and the two
# formats. The sources are index.md here, ROADMAP.md, and
# the files of the two format folders: SPEC.md and RATIONALE.md, rendered as
# they are, and overview.qmd, which is run. The folders are copied to a work
# folder and rendered there as one Quarto website, so nothing is written into
# them and no HTML of the site is kept in the repository.
#
# The overviews report the results of the verification suites, so run
# verify/run.sh in both folders first. Rendering needs Quarto, R with knitr and
# reticulate, and the Python named in SCEPTRE_FORMAT_PYTHON (default: python3). The analysis
# format's overview reads an example into a sceptre_object, so R also needs
# sceptre 0.99.1. Run from anywhere.
set -e
cd "$(dirname "$0")"
export LANG=en_US.UTF-8
export RETICULATE_PYTHON="${SCEPTRE_FORMAT_PYTHON:-$(command -v python3)}"
Rscript -e 'if (!requireNamespace("sceptre", quietly = TRUE) || packageVersion("sceptre") != "0.99.1") stop("The overview of the analysis format needs sceptre 0.99.1.", call. = FALSE)'
rm -rf _work _site
mkdir _work
cp _quarto.yml index.md rules.lua site.css _work/
# the roadmap: its heading becomes the page title, its diagram is run, and its
# links to the format folders go to their overviews
[ "$(head -n 1 ../ROADMAP.md)" = "# SCEPTRE roadmap" ] || { echo "ROADMAP.md does not start with its heading." >&2; exit 1; }
printf -- '---\ntitle: Roadmap\n---\n' > _work/roadmap.qmd
sed -e '1d' -e 's/^```mermaid$/```{mermaid}/' \
    -e 's#(perturb-seq-data-format/)#(perturb-seq-data-format/overview.qmd)#g' \
    -e 's#(sceptre-analysis-format/)#(sceptre-analysis-format/overview.qmd)#g' ../ROADMAP.md >> _work/roadmap.qmd
for format in perturb-seq-data-format sceptre-analysis-format; do
    if [ ! -d "../$format/verify/out/results" ]; then
        echo "Run $format/verify/run.sh first: the overview reports its results." >&2
        exit 1
    fi
    mkdir "_work/$format"
    (cd "../$format" && tar -cf - --exclude ./rlib .) | (cd "_work/$format" && tar -xf -)
    [ -d "../$format/rlib" ] && ln -s "$(cd "../$format/rlib" && pwd)" "_work/$format/rlib"
done
# the home page's badges carry each format's version, as its specification states it
for format in perturb-seq-data-format sceptre-analysis-format; do
    version="$(sed -n 's/^\*\*Version \(.*\)\.\*\*$/\1/p' "../$format/SPEC.md" | head -n 1)"
    [ -n "$version" ] || { echo "No version line in $format/SPEC.md." >&2; exit 1; }
    label="$(printf 'version %s' "${version%%,*}" | sed -e 's/ /%20/g' -e 's/(/%28/g' -e 's/)/%29/g')"
    sed -i.bak "s/{{version $format}}/$label/" _work/index.md
done
rm _work/index.md.bak
quarto render _work
mv _work/_site _site
rm -rf _work
