# Libraries for the format code: nothing but standard packages. `nanoparquet`
# and `jsonvalidate` are taken from ./rlib if they are there (see README.md).
# Run from the folder that holds format/ and spec/.
format_dir <- Sys.getenv("SCEPTRE_FORMAT_DIR", getwd())
if (dir.exists(file.path(format_dir, "rlib"))) .libPaths(c(file.path(format_dir, "rlib"), .libPaths()))
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(Matrix)
})
source(file.path(format_dir, "format", "sceptre_format.R"))
