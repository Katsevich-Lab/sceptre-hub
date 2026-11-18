# Libraries for the format code: nothing but standard packages. `nanoparquet`
# and `jsonvalidate` are taken from ./rlib if they are there (see README.md).
# Run from the folder that holds format/ and spec/.
proposal_dir <- Sys.getenv("SCEPTRE_FORMAT_DIR", getwd())
if (dir.exists(file.path(proposal_dir, "rlib"))) .libPaths(c(file.path(proposal_dir, "rlib"), .libPaths()))
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(Matrix)
})
source(file.path(proposal_dir, "format", "sceptre_format.R"))
