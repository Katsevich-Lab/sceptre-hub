# Libraries and code for the analysis format: the data file format's code
# first, then this folder's. Run from this folder's parent (the one that holds
# format/ and spec/), or name it in SCEPTRE_ANALYSIS_FORMAT_DIR. The data file
# format is taken from ../perturb-seq-data-format, or from
# PERTURB_SEQ_DATA_FORMAT_DIR.
analysis_format_dir <- Sys.getenv("SCEPTRE_ANALYSIS_FORMAT_DIR", getwd())
Sys.setenv(SCEPTRE_ANALYSIS_FORMAT_DIR = analysis_format_dir)
if (!nzchar(Sys.getenv("PERTURB_SEQ_DATA_FORMAT_DIR"))) {
    Sys.setenv(PERTURB_SEQ_DATA_FORMAT_DIR = normalizePath(file.path(analysis_format_dir, "..", "perturb-seq-data-format")))
}
source(file.path(Sys.getenv("PERTURB_SEQ_DATA_FORMAT_DIR"), "format", "setup.R"))
source(file.path(analysis_format_dir, "format", "sceptre_analysis.R"))
