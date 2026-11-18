# Libraries for the prototype: nanoparquet from ./rlib, and the sceptre3 build of
# sceptre from the tour library (see README.md). Run from the proposal folder.
proposal_dir <- Sys.getenv("SCEPTRE_FORMAT_DIR", getwd())
.libPaths(c(
    file.path(proposal_dir, "rlib"),
    Sys.getenv("SCEPTRE_TOUR_RLIB", path.expand("~/code/research/sceptre-tour-rlib")),
    .libPaths()
))
suppressPackageStartupMessages({
    library(sceptre)
    library(SingleCellExperiment)
})
source(file.path(proposal_dir, "format", "sceptre_format.R"))
