# The R half of analysis_scale.py: one measurement per call, answer as JSON.
#   Rscript analysis_scale.R read_all <nanoparquet|arrow|arrow_table> <file>
#   Rscript analysis_scale.R read_where arrow <file> <column> <value> [<response_id>]
#   Rscript analysis_scale.R fits_read_each nanoparquet <folder>
# `arrow` gives a data frame; `arrow_table` leaves the data in Arrow's own memory.
arguments <- commandArgs(trailingOnly = TRUE)
operation <- arguments[1]; reader <- arguments[2]; path <- arguments[3]
start <- proc.time()[["elapsed"]]
if (operation == "read_all") {
  table <- switch(reader,
    nanoparquet = nanoparquet::read_parquet(path),
    arrow = arrow::read_parquet(path),
    arrow_table = arrow::read_parquet(path, as_data_frame = FALSE))
  rows <- nrow(table)
} else if (operation == "read_where") {
  suppressPackageStartupMessages(library(dplyr))
  rows <- arrow::open_dataset(path)
  rows <- if (arguments[4] == "grna_target") filter(rows, grna_target == arguments[5]) else filter(rows, response_id == arguments[5])
  if (length(arguments) == 6) rows <- filter(rows, response_id == arguments[6])
  rows <- nrow(collect(rows))
} else {
  files <- sort(list.files(path, pattern = "\\.parquet$", full.names = TRUE))
  for (f in files) {
    meta <- nanoparquet::read_parquet_metadata(f)$file_meta_data$key_value_metadata[[1]]
    stopifnot("sceptre_fit_identity" %in% meta$key)
  }
  rows <- nrow(do.call(rbind, lapply(files, nanoparquet::read_parquet)))
}
cat(sprintf('{"rows": %d, "seconds": %.2f}\n', rows, proc.time()[["elapsed"]] - start))
