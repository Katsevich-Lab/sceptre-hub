# Libraries for the data file format's code: nothing but standard packages.
# `nanoparquet` and `jsonvalidate` are taken from ./rlib if they are there (see
# README.md). Run from this folder's parent (the one that holds format/ and
# spec/), or name it in PERTURB_SEQ_DATA_FORMAT_DIR.
data_format_dir <- Sys.getenv("PERTURB_SEQ_DATA_FORMAT_DIR", getwd())
Sys.setenv(PERTURB_SEQ_DATA_FORMAT_DIR = data_format_dir)
if (dir.exists(file.path(data_format_dir, "rlib"))) .libPaths(c(file.path(data_format_dir, "rlib"), .libPaths()))
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(Matrix)
})
source(file.path(data_format_dir, "format", "perturb_seq_data.R"))
