# The sceptre file format (draft, format 0.1): R reader and writer.
#
# SPEC.md states the rules and RATIONALE.md the reasons. Section numbers below
# (S4, T3, ...) refer to both.
#
# Two layers:
#   1. Files <-> standard R containers, with no use of the sceptre package:
#        fmt_write_screen(), fmt_read_screen()   data file <-> SingleCellExperiment
#        fmt_init_analysis(), fmt_write_stage(), fmt_read_analysis(), fmt_clear_stage()
#   2. Files <-> sceptre_object (needs the sceptre3 build of sceptre):
#        write_sceptre_screen(), write_sceptre_analysis(), read_sceptre_object()
#   3. Checking files against SPEC.md, with no use of the sceptre package and
#      without building any object:
#        validate_screen(), validate_analysis()
#
# The tables of an analysis, the stages and the shape of the JSON records are
# not written here: they are read from spec/ (tables.json, screen.json and two
# JSON Schemas), which the Python library reads too.
#
# Names. The data file names a modality for what was measured (`rna`), and the
# gRNA modality is always `grna` (S2). An analysis says which modality is its
# response (A1). The sceptre package calls the response modality `response`
# whatever was measured, so layer 2 translates covariate names, formulas and
# the gRNA table's coordinate columns between the two.

# The format version, "MAJOR.MINOR" (S10.3). A minor increase adds things an
# older reader can ignore; a major increase is one an older reader must refuse.
FORMAT_VERSION <- "0.1"
format_parts <- function(version) {
    version <- as.character(as.vector(version))
    if (length(version) != 1L || is.na(version) || !grepl("^[0-9]+\\.[0-9]+$", version)) return(NULL)
    as.integer(strsplit(version, ".", fixed = TRUE)[[1]])
}
known_major <- function(version) identical(format_parts(version)[1], format_parts(FORMAT_VERSION)[1])
SPEC_DIR <- file.path(Sys.getenv("SCEPTRE_FORMAT_DIR", getwd()), "spec")
TABLE_SPEC <- jsonlite::fromJSON(file.path(SPEC_DIR, "tables.json"), simplifyVector = FALSE)
SCREEN_SPEC <- jsonlite::fromJSON(file.path(SPEC_DIR, "screen.json"), simplifyVector = FALSE)
GRNA <- "grna"
DEFAULT_RESPONSE_MODALITY <- "rna"
# the gRNA table's target coordinates: name in the package = name in the file (S8)
GRNA_COORDINATES <- c(chr = "target_chr", start = "target_start", end = "target_end")
STAGES <- unlist(TABLE_SPEC$stages)
# calibration, power and discovery do not depend on one another (A4)
STAGE_RANK <- vapply(TABLE_SPEC$stage_rank, as.integer, integer(1))[STAGES]
STAGE_FUNCTIONS <- c(
    parameters = "set_analysis_parameters", assignment = "assign_grnas", qc = "run_qc",
    calibration = "run_calibration_check", power = "run_power_check", discovery = "run_discovery_analysis"
)
# The shared list of response-fitting methods (T6). A method fitted in one
# language only is never reused by the other.
FIT_METHODS <- list(batched = c("R", "Python"), sceptre = c("R", "Python"), pysceptre = c("R", "Python"), glmGamPoi = "R")
ROW_GROUP_ROWS <- as.integer(TABLE_SPEC$max_rows_per_row_group)

fmt_stop <- function(...) stop(..., call. = FALSE)

# ---------------------------------------------------------------- small helpers

sha256_ids <- function(ids) unname(tools::sha256sum(bytes = charToRaw(paste(enc2utf8(ids), collapse = "\n"))))
new_id <- function() unname(tools::md5sum(bytes = charToRaw(paste(Sys.time(), Sys.getpid(), tempfile(), stats::runif(1)))))
byte_order <- function(table, key) do.call(order, c(unname(as.list(table[key])), list(method = "radix")))

# Each file is written under a temporary name and renamed into place (A6).
atomically <- function(path, write) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- paste0(path, ".tmp")
    write(temporary)
    if (!file.rename(temporary, path)) fmt_stop("Could not write ", path)
    invisible(path)
}

write_json <- function(value, path) {
    atomically(path, function(temporary) {
        # 17 significant digits: a double read back is the double written
        writeLines(jsonlite::toJSON(value, auto_unbox = TRUE, digits = I(17), pretty = TRUE, null = "null", na = "null"), temporary)
    })
}
read_json <- function(path) {
    jsonlite::fromJSON(paste(readLines(path, warn = FALSE), collapse = "\n"), simplifyVector = TRUE, simplifyDataFrame = FALSE)
}

# ---------------------------------------------------------------- tables (T0)

# Rows sorted by the key in byte order; row groups of at most a million rows.
write_table <- function(table, path, key) {
    table <- as.data.frame(table, stringsAsFactors = FALSE)
    factors <- vapply(table, is.factor, logical(1))
    table[factors] <- lapply(table[factors], as.character)
    table <- table[byte_order(table, key), , drop = FALSE]
    rownames(table) <- NULL
    atomically(path, function(temporary) {
        nanoparquet::write_parquet(
            table, temporary,
            options = nanoparquet::parquet_options(num_rows_per_row_group = ROW_GROUP_ROWS)
        )
    })
    nrow(table)
}

table_files <- function(path) {
    if (dir.exists(path)) sort(list.files(path, pattern = "\\.parquet$", full.names = TRUE)) else path
}
table_rows <- function(path) sum(vapply(table_files(path), function(f) nanoparquet::read_parquet_info(f)$num_rows, numeric(1)))

# One .parquet file or a folder of part files; only the known columns where
# `columns` is given; identifiers as text whether stored plain or as a dictionary.
read_table <- function(path, columns = NULL) {
    table <- do.call(rbind, lapply(table_files(path), function(file) {
        part <- as.data.frame(nanoparquet::read_parquet(file), stringsAsFactors = FALSE)
        class(part) <- "data.frame"
        if (!is.null(columns)) part <- part[intersect(columns, names(part))]
        factors <- vapply(part, is.factor, logical(1))
        part[factors] <- lapply(part[factors], as.character)
        part
    }))
    rownames(table) <- NULL
    table
}

# Every table of an analysis, from spec/tables.json: the stage that writes it,
# its key, and the columns a reader takes (NULL: all of them).
unit_key <- function(table) intersect(c("grna_target", "grna_id", "response_id"), names(table))
TABLES <- stats::setNames(lapply(TABLE_SPEC$tables, function(definition) {
    columns <- vapply(definition$columns, function(column) column$name, character(1))
    list(
        stage = definition$stage, key = unlist(definition$key), required = definition$required,
        columns = if (identical(definition$extra_columns, "ignored")) columns else NULL,
        definition = definition
    )
}), vapply(TABLE_SPEC$tables, function(definition) definition$name, character(1)))
# a key column that is not required is part of the key when it is present
table_key <- function(name, table) intersect(TABLES[[name]]$key, names(table))

# ---------------------------------------------------------------- the data file

# Column types of the cell and feature tables (S7).
#   covariates (cell tables): text and true/false become categoricals
#   feature tables:           text stays text, as a categorical only where it has
#                             missing values (R's MuData cannot store a missing string)
#   both:                     an integer column with missing values becomes doubles
as_file_frame <- function(frame, covariates) {
    frame <- as.data.frame(frame, optional = TRUE, stringsAsFactors = FALSE)
    for (name in names(frame)) {
        column <- frame[[name]]
        if (is.logical(column)) {
            frame[[name]] <- factor(column, levels = c(FALSE, TRUE), labels = c("FALSE", "TRUE"))
        } else if (is.character(column) && (covariates || anyNA(column))) {
            frame[[name]] <- factor(column)
        } else if (is.integer(column) && anyNA(column)) {
            frame[[name]] <- as.double(column)
        }
    }
    frame
}
from_file_frame <- function(frame, covariates) {
    frame <- as.data.frame(frame, optional = TRUE, stringsAsFactors = FALSE)
    for (name in names(frame)) {
        column <- frame[[name]]
        if (is.factor(column) && !covariates) frame[[name]] <- as.character(column)
        if (is.double(column)) frame[[name]] <- replace(column, is.nan(column), NA)
    }
    frame
}

# The model's covariate table is split in the file (S5): facts about the cell
# are in the shared cell table, measurements of a modality in that modality's,
# without the prefix `<modality>_`.
split_covariates <- function(covariates, modalities) {
    covariates <- as.data.frame(covariates, optional = TRUE, stringsAsFactors = FALSE)
    of <- function(modality) startsWith(names(covariates), paste0(modality, "_"))
    out <- lapply(stats::setNames(modalities, modalities), function(modality) {
        part <- covariates[of(modality)]
        names(part) <- substring(names(part), nchar(modality) + 2L)
        part
    })
    out$shared <- covariates[!Reduce(`|`, lapply(modalities, of))]
    out
}
join_covariates <- function(parts, modalities) {
    prefixed <- lapply(modalities, function(modality) {
        part <- parts[[modality]]
        if (ncol(part) > 0L) names(part) <- paste0(modality, "_", names(part))
        part
    })
    frames <- c(prefixed, list(parts$shared))
    frames <- frames[vapply(frames, ncol, integer(1)) > 0L]
    if (length(frames) == 0L) parts$shared else do.call(cbind, frames)
}

scalar <- function(x) {
    attr(x, "encoding-scalar") <- TRUE
    x
}

# Counts as 32-bit integers, cell positions as 32-bit integers, column pointers
# as 64-bit integers, uncompressed (S4, S11). MuData writes doubles and 32-bit
# pointers; this rewrites the three arrays of one modality.
write_count_types <- function(file, modality) {
    group <- paste0("mod/", modality, "/X/")
    replace <- function(name, h5type) {
        values <- as.vector(rhdf5::h5read(file, paste0(group, name), bit64conversion = "double"))
        if (name == "data" && any(values != round(values))) fmt_stop("Counts must be integers.")
        rhdf5::h5delete(file, paste0(group, name))
        rhdf5::h5createDataset(file, paste0(group, name), dims = length(values), H5type = h5type, chunk = NULL, level = 0)
        rhdf5::h5write(if (h5type == "H5T_STD_I64LE") as.numeric(values) else as.integer(values), file, paste0(group, name))
    }
    replace("data", "H5T_STD_I32LE")
    replace("indices", "H5T_STD_I32LE")
    replace("indptr", "H5T_STD_I64LE")
}

# The modality an analysis tests, among those of a data file: the one asked for,
# or the only one besides `grna` (S2).
response_modality_of <- function(modalities, requested = NULL, what = "The data file") {
    if (!GRNA %in% modalities) fmt_stop(what, " has no modality named '", GRNA, "'.")
    candidates <- setdiff(modalities, GRNA)
    if (is.null(requested)) {
        if (length(candidates) == 1L) return(candidates)
        fmt_stop(what, " has the modalities ", paste(candidates, collapse = ", "), " besides '", GRNA, "'; say which one is the response.")
    }
    if (!requested %in% candidates) fmt_stop(what, " has no modality named '", requested, "'.")
    requested
}

# Write a SingleCellExperiment as a data file.
#   assay "counts"          counts of the measured modality, features x cells
#   colData                 the model's covariate table, a modality's columns
#                           prefixed with its name (`rna_n_umis`, `grna_n_umis`)
#   rowData                 feature metadata
#   altExp "gRNA"           gRNA counts; its rowData is the gRNA table (grna_target, ...)
#   other altExps           further modalities over the same cells, written under
#                           their names; analyses that do not name them ignore them
#   metadata(sce)$sceptre   list(moi = "high" or "low", response_modality = the
#                           main experiment's modality name, "rna" if absent);
#                           `id` is kept if present
fmt_write_screen <- function(sce, file) {
    if (file.exists(file)) fmt_stop(file, " exists; a data file is written once.")
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    stamp <- S4Vectors::metadata(sce)$sceptre
    if (!isTRUE(stamp$moi %in% c("high", "low"))) fmt_stop("metadata(sce)$sceptre$moi must be 'high' or 'low'.")
    main <- if (is.null(stamp$response_modality)) DEFAULT_RESPONSE_MODALITY else stamp$response_modality
    others <- setdiff(SingleCellExperiment::altExpNames(sce), "gRNA")
    names_used <- c(main, GRNA, others)
    if (anyDuplicated(names_used) > 0L || !all(grepl("^[A-Za-z][A-Za-z0-9_]*$", names_used))) {
        fmt_stop("Modality names must be distinct, start with a letter and use letters, digits and underscores; '", GRNA, "' is the gRNA modality.")
    }
    grna <- SingleCellExperiment::altExp(sce, "gRNA")
    # cells without identifiers get their positions (S9)
    cells <- if (is.null(colnames(sce))) as.character(seq_len(ncol(sce))) else colnames(sce)
    features <- c(rownames(sce), rownames(grna), unlist(lapply(others, function(name) rownames(SingleCellExperiment::altExp(sce, name)))))
    if (anyDuplicated(cells) > 0L) fmt_stop("Cell identifiers must be unique.")
    if (anyDuplicated(features) > 0L) fmt_stop("Feature identifiers must be unique across the modalities of a file.")
    grna_table <- SummarizedExperiment::rowData(grna)
    if (is.null(grna_table$grna_target) || anyNA(grna_table$grna_target)) fmt_stop("Every gRNA needs a `grna_target`.")
    covariates <- split_covariates(SummarizedExperiment::colData(sce), c(main, GRNA))
    modality <- function(experiment, obs, var) {
        # every modality's cell table has its total count and its number of
        # features with a nonzero count; a writer computes them if they are absent (S6)
        counts <- SummarizedExperiment::assay(experiment, "counts")
        obs <- as.data.frame(obs, optional = TRUE, stringsAsFactors = FALSE)
        if (is.null(obs$n_umis)) obs$n_umis <- as.double(Matrix::colSums(counts))
        if (is.null(obs$n_nonzero)) obs$n_nonzero <- as.integer(Matrix::colSums(counts != 0))
        SingleCellExperiment::SingleCellExperiment(
            # row-compressed (features x cells) in R is column-compressed
            # (cells x features) in the file: one feature is one contiguous read (S3)
            assays = list(counts = methods::as(
                methods::as(SummarizedExperiment::assay(experiment, "counts"), "generalMatrix"), "RsparseMatrix"
            )),
            colData = S4Vectors::DataFrame(as_file_frame(obs, TRUE), check.names = FALSE, row.names = cells),
            rowData = S4Vectors::DataFrame(as_file_frame(var, FALSE), check.names = FALSE, row.names = rownames(experiment))
        )
    }
    written <- FALSE
    on.exit(if (!written) unlink(file), add = TRUE)
    experiments <- c(
        stats::setNames(list(
            modality(sce, covariates[[main]], SummarizedExperiment::rowData(sce)),
            modality(grna, covariates[[GRNA]], grna_table)
        ), c(main, GRNA)),
        lapply(stats::setNames(others, others), function(name) {
            other <- SingleCellExperiment::altExp(sce, name)
            modality(other, data.frame(row.names = colnames(other)), SummarizedExperiment::rowData(other))
        })
    )
    MuData::writeH5MU(
        MultiAssayExperiment::MultiAssayExperiment(
            experiments,
            colData = S4Vectors::DataFrame(as_file_frame(covariates$shared, TRUE), check.names = FALSE, row.names = cells),
            # a native record with scalar entries (S10)
            metadata = list(sceptre = list(
                format = scalar(FORMAT_VERSION), moi = scalar(stamp$moi),
                id = scalar(if (is.null(stamp$id)) new_id() else stamp$id)
            ))
        ),
        file
    )
    for (name in names(experiments)) write_count_types(file, name)
    written <- TRUE
    invisible(file)
}

check_format <- function(format, what) {
    if (is.null(format_parts(format))) fmt_stop(what, " has no sceptre format version.")
    if (!known_major(format)) {
        fmt_stop(what, " is in format ", format, "; this reader knows format ", FORMAT_VERSION, " and its minor versions.")
    }
}

# Read a data file as a SingleCellExperiment (the layout fmt_write_screen takes).
# `response_modality` picks the modality that becomes the main experiment; it
# may be left out when the file has one modality besides `grna`. Modalities the
# reader is not asked for are ignored.
fmt_read_screen <- function(file, response_modality = NULL) {
    mae <- withCallingHandlers(
        MuData::readH5MU(file),
        warning = function(w) {
            if (grepl("sampleMap", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
        }
    )
    stamp <- lapply(S4Vectors::metadata(mae)$sceptre, function(x) as.vector(x)[1])
    check_format(stamp$format, file)
    main <- response_modality_of(names(mae), response_modality, file)
    cells <- colnames(mae[[main]])
    if (!identical(cells, colnames(mae[[GRNA]]))) fmt_stop("The modalities must hold the same cells in the same order.")
    cell_table <- function(frame) {
        frame <- from_file_frame(frame, TRUE)
        frame <- frame[match(cells, rownames(frame)), , drop = FALSE]
        rownames(frame) <- cells
        frame
    }
    covariates <- join_covariates(c(
        stats::setNames(list(
            cell_table(SummarizedExperiment::colData(mae[[main]])),
            cell_table(SummarizedExperiment::colData(mae[[GRNA]]))
        ), c(main, GRNA)),
        list(shared = cell_table(SummarizedExperiment::colData(mae)))
    ), c(main, GRNA))
    experiment <- function(name) {
        e <- mae[[name]]
        SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = methods::as(SummarizedExperiment::assay(e, 1L), "CsparseMatrix")),
            rowData = S4Vectors::DataFrame(
                from_file_frame(SummarizedExperiment::rowData(e), FALSE), check.names = FALSE, row.names = rownames(e)
            )
        )
    }
    sce <- experiment(main)
    SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(covariates, check.names = FALSE, row.names = cells)
    SingleCellExperiment::altExp(sce, "gRNA") <- experiment(GRNA)
    S4Vectors::metadata(sce) <- list(sceptre = list(
        format = as.character(stamp$format), moi = stamp$moi, id = stamp$id, response_modality = main
    ))
    sce
}

# What an analysis records about its data file (A1).
screen_link <- function(sce, file = NULL, directory = NULL) {
    grna <- SingleCellExperiment::altExp(sce, "gRNA")
    list(
        id = S4Vectors::metadata(sce)$sceptre$id,
        n_cells = ncol(sce), n_responses = nrow(sce), n_grnas = nrow(grna),
        cell_ids_sha256 = sha256_ids(colnames(sce)),
        response_ids_sha256 = sha256_ids(rownames(sce)),
        grna_ids_sha256 = sha256_ids(rownames(grna)),
        # a hint only; the identity above is what a reader checks
        path = if (is.null(file)) NULL else relative_path(file, directory)
    )
}
relative_path <- function(file, directory) {
    file <- strsplit(normalizePath(file, mustWork = FALSE), "/", fixed = TRUE)[[1]]
    from <- strsplit(normalizePath(directory, mustWork = FALSE), "/", fixed = TRUE)[[1]]
    shared <- sum(cumprod(utils::head(file, min(length(file), length(from))) == utils::head(from, min(length(file), length(from)))))
    paste(c(rep("..", length(from) - shared), file[-seq_len(shared)]), collapse = "/")
}

# ---------------------------------------------------------------- the analysis folder

# Start an analysis folder: which data file it analyses and which of its
# modalities is the response (A1).
fmt_init_analysis <- function(directory, link, response_modality) {
    if (length(list.files(directory)) > 0L) fmt_stop(directory, " is not empty.")
    write_json(
        list(format = FORMAT_VERSION, response_modality = response_modality, screen = link),
        file.path(directory, "analysis.json")
    )
    invisible(directory)
}
# The modality an analysis folder names as its response.
fmt_response_modality <- function(directory) {
    name <- read_json(file.path(directory, "analysis.json"))$response_modality
    if (is.null(name)) fmt_stop(directory, " does not say which modality is its response.")
    name
}

stage_file <- function(directory, stage) file.path(directory, "stages", paste0(stage, ".json"))
fmt_stages_run <- function(directory) STAGES[file.exists(stage_file(directory, STAGES))]

# Write one stage: its tables first, then its record, which lists them (A2, A5).
# `record` holds settings, and optionally summaries and anything else of the stage.
fmt_write_stage <- function(directory, stage, record, tables = list(), provenance = NULL) {
    if (!file.exists(file.path(directory, "analysis.json"))) fmt_stop(directory, " is not an analysis folder.")
    earlier <- names(STAGE_RANK)[STAGE_RANK < STAGE_RANK[[stage]]]
    missing <- setdiff(earlier, fmt_stages_run(directory))
    if (length(missing) > 0L) fmt_stop("Stage '", stage, "' needs the stage(s) ", paste(missing, collapse = ", "), " first.")
    if (file.exists(stage_file(directory, stage))) fmt_stop("Stage '", stage, "' is already written; clear it first.")
    files <- lapply(names(tables), function(name) {
        if (is.null(TABLES[[name]]) || TABLES[[name]]$stage != stage) fmt_stop("Stage '", stage, "' does not write the table ", name, ".")
        path <- paste0(name, ".parquet")
        list(path = path, rows = write_table(tables[[name]], file.path(directory, path), table_key(name, tables[[name]])))
    })
    if (is.null(provenance)) provenance <- list(package = "sceptre_format.R", version = "prototype", language = "R")
    if (is.null(provenance$written_at)) provenance$written_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    write_json(c(list(stage = stage), record, list(provenance = provenance, files = files)), stage_file(directory, stage))
    invisible(directory)
}

# Remove a stage and every later stage: their records and their files (A4).
fmt_clear_stage <- function(directory, stage) {
    for (later in names(STAGE_RANK)[STAGE_RANK > STAGE_RANK[[stage]] | names(STAGE_RANK) == stage]) {
        path <- stage_file(directory, later)
        if (!file.exists(path)) next
        record <- read_json(path)
        # the record goes first: a stage without a record has not run
        unlink(path)
        unlink(file.path(directory, vapply(record$files, function(f) f$path, character(1))), recursive = TRUE)
    }
    invisible(directory)
}

# Read an analysis folder: its link, the record of each stage that ran, and
# their tables (named as in TABLES).
fmt_read_analysis <- function(directory) {
    analysis <- read_json(file.path(directory, "analysis.json"))
    check_format(analysis$format, directory)
    stages <- lapply(stats::setNames(nm = fmt_stages_run(directory)), function(stage) read_json(stage_file(directory, stage)))
    tables <- list()
    for (record in stages) {
        for (file in record$files) {
            path <- file.path(directory, file$path)
            if (!file.exists(path)) fmt_stop("The stage '", record$stage, "' lists ", file$path, ", which is missing.")
            if (table_rows(path) != file$rows) {
                fmt_stop(file$path, " has ", table_rows(path), " rows; the stage '", record$stage, "' wrote ", file$rows, ".")
            }
            name <- sub("\\.parquet$", "", file$path)
            tables[[name]] <- read_table(path, TABLES[[name]]$columns)
        }
    }
    list(
        format = analysis$format, response_modality = analysis$response_modality,
        screen = analysis$screen, stages = stages, tables = tables
    )
}

# Refuse a data file other than the one the analysis was run on (A1).
check_link <- function(recorded, sce, what) {
    current <- screen_link(sce)
    for (field in c("id", "cell_ids_sha256", "response_ids_sha256", "grna_ids_sha256")) {
        if (!identical(recorded[[field]], current[[field]])) {
            fmt_stop("The analysis ", what, " was run on another data file (its `", field, "` differs).")
        }
    }
}

# ---------------------------------------------------------------- derivations

# Values a package holds that follow from the files; both languages derive them
# the same way, and the sceptre_object writer checks each against the object.

# Each cell's gRNA with the most UMIs (the first on a tie) and its share of the
# cell's gRNA UMIs (S6).
top_grna <- function(grna_counts) {
    counts <- methods::as(methods::as(grna_counts, "generalMatrix"), "TsparseMatrix")
    entries <- data.frame(grna = counts@i + 1L, cell = counts@j + 1L, umis = counts@x)
    entries <- entries[order(entries$cell, -entries$umis, entries$grna), ]
    top <- entries[!duplicated(entries$cell), ]
    max_grna <- rep(NA_character_, ncol(grna_counts))
    max_grna[top$cell] <- rownames(grna_counts)[top$grna]
    fraction <- rep(NA_real_, ncol(grna_counts))
    fraction[top$cell] <- top$umis / Matrix::colSums(grna_counts)[top$cell]
    list(max_grna = max_grna, max_grna_frac_umis = fraction)
}

# The gRNA group of each gRNA under an integration strategy.
grna_groups <- function(grna_table, strategy) {
    if (identical(strategy, "union")) return(grna_table$grna_target)
    ifelse(grna_table$grna_target == "non-targeting", "non-targeting", grna_table$grna_id)
}

# The pairs as a package tests them: one row per (gRNA group, response).
expand_pairs <- function(pairs, grna_table) {
    groups <- unique(grna_table[grna_table$grna_group != "non-targeting", c("grna_target", "grna_group")])
    out <- merge(pairs, groups, by = "grna_target", sort = FALSE)
    out[c("grna_target", "response_id", "grna_group")]
}

# Each gRNA's cells (positions among the data file's cells, ascending), from
# the assignment table (T1).
cells_by_grna <- function(assignments, grna_ids, cells) {
    position <- match(assignments$cell_id, cells)
    by_grna <- split(position, factor(assignments$grna_id, levels = grna_ids))
    lapply(by_grna, sort)
}
# Each targeting group's cells (gRNA by gRNA, each cell once) and each
# non-targeting gRNA's cells.
group_cells <- function(initial, grna_table) {
    targeting <- grna_table$grna_group != "non-targeting"
    groups <- unique(grna_table$grna_group[targeting])
    list(
        grna_group_idxs = lapply(stats::setNames(groups, groups), function(group) {
            unique(unlist(initial[grna_table$grna_id[grna_table$grna_group == group]], use.names = FALSE))
        }),
        indiv_nt_grna_idxs = initial[grna_table$grna_id[!targeting]]
    )
}
grnas_per_cell <- function(initial, n_cells) tabulate(unlist(initial, use.names = FALSE), nbins = n_cells)
# In low MOI, the cells that do not carry exactly one gRNA.
not_one_grna <- function(initial, n_cells) which(grnas_per_cell(initial, n_cells) != 1L)

# After QC: positions among the cells in use; with the NT cells as the control
# group, each non-targeting gRNA's cells as positions among all NT cells.
qc_assignments <- function(raw, cells_in_use, nt_cells_control) {
    among <- function(cells) {
        position <- match(cells, cells_in_use)
        sort(position[!is.na(position)])
    }
    out <- lapply(raw, function(group) lapply(group, among))
    if (nt_cells_control) {
        out$all_nt_idxs <- unique(unlist(out$indiv_nt_grna_idxs, use.names = FALSE))
        stop <- cumsum(lengths(out$indiv_nt_grna_idxs))
        start <- c(0L, utils::head(stop, -1L)) + 1L
        out$indiv_nt_grna_idxs <- stats::setNames(
            lapply(seq_along(stop), function(i) seq(start[i], stop[i])), names(out$indiv_nt_grna_idxs)
        )
    }
    out
}

# How many pairs pass pairwise QC, as sceptre counts them: under the bonferroni
# strategy a (target, response) pair passes if any of its gRNAs does.
n_pairs_passing <- function(pair_qc, strategy) {
    if (!identical(strategy, "bonferroni")) return(sum(pair_qc$pass_qc))
    sum(tapply(pair_qc$pass_qc, paste(pair_qc$grna_target, pair_qc$response_id, sep = "\r"), any))
}

# The largest p-value called significant, or NULL: `significant` is then
# p_value <= threshold, which no part file needs the whole table for (T5).
significance_threshold <- function(result) {
    called <- which(result$significant)
    if (length(called) == 0L) NULL else max(result$p_value[called])
}
add_significant <- function(result, threshold) {
    result$significant <- if (is.null(threshold)) ifelse(is.na(result$p_value), NA, FALSE) else result$p_value <= threshold
    result
}

# `log_2_fold_change` follows from `fold_change`, so a results table that has
# `fold_change` does not store it (T5). The formula is sceptre's, to the last
# bit: log(x) / log(2), which is not always what log2(x) returns.
log_2_of <- function(fold_change) log(fold_change) / log(2)
log_2_follows <- function(result) {
    if (is.null(result$fold_change) || is.null(result$log_2_fold_change)) return(FALSE)
    derived <- log_2_of(result$fold_change)
    held <- result$log_2_fold_change
    all((is.na(derived) & is.na(held)) | (!is.na(derived) & !is.na(held) & derived == held))
}
add_log_2_fold_change <- function(result) {
    if (!is.null(result$log_2_fold_change) || is.null(result$fold_change)) return(result)
    after <- max(match(c("fold_change", "se_fold_change"), names(result)), na.rm = TRUE)
    columns <- append(as.list(result), list(log_2_fold_change = log_2_of(result$fold_change)), after = after)
    as.data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE)
}

# Resampled statistics: packages hold them as columns z_null_1, z_null_2, ...
# padded with missing values; the file holds one list per pair (T5).
resample_columns <- function(table) grep("^z_null_[0-9]+$", names(table), value = TRUE)
resamples_table <- function(result) {
    columns <- resample_columns(result)
    if (length(columns) == 0L) return(NULL)
    values <- as.matrix(result[columns])
    out <- result[unit_key(result)]
    out$z_null <- lapply(seq_len(nrow(values)), function(i) unname(values[i, !is.na(values[i, ])]))
    out
}
add_resamples <- function(result, resamples) {
    if (is.null(resamples)) return(result)
    key <- function(table) do.call(paste, c(unname(as.list(table[unit_key(table)])), sep = "\r"))
    lists <- resamples$z_null[match(key(result), key(resamples))]
    width <- max(lengths(lists))
    values <- t(vapply(lists, function(z) c(z, rep(NA_real_, width - length(z))), numeric(width)))
    colnames(values) <- paste0("z_null_", seq_len(width))
    cbind(result, as.data.frame(values))
}

# Additive formulas only: each term a covariate name, log(x) or log(x + 1),
# with 0 or 1 for the intercept (T7).
formula_in_subset <- function(text) {
    rhs <- gsub("[[:space:]]", "", sub("^[^~]*~", "", text))
    characters <- strsplit(rhs, "")[[1]]
    depth <- cumsum(characters == "(") - cumsum(characters == ")")
    plus <- which(characters == "+" & depth == 0L)
    terms <- substring(rhs, c(1L, plus + 1L), c(plus - 1L, nchar(rhs)))
    name <- "[A-Za-z][A-Za-z0-9_.]*"
    all(grepl(paste0("^(0|1|", name, "|log\\(", name, "\\)|log\\(", name, "\\+1\\))$"), terms))
}
# SHA-256 of the matrix's values as little-endian doubles, column by column.
design_sha256 <- function(matrix) {
    unname(tools::sha256sum(bytes = writeBin(as.double(matrix), raw(), endian = "little")))
}

# ---------------------------------------------------------------- sceptre_object

sceptre_ns <- function() asNamespace("sceptre")
field <- function(object, name) sceptre_ns()$object_field(object, name)
`field<-` <- function(object, name, value) {
    setter <- get("object_field<-", envir = sceptre_ns())
    setter(object, name, value)
}

# The package's names <-> the file's names. The package prefixes the response
# modality's covariates with `response_`; the file uses the modality's name.
swap_prefix <- function(text, from, to) {
    gsub(paste0("(?<![A-Za-z0-9_.])", from, "_"), paste0(to, "_"), text, perl = TRUE)
}
to_file_names <- function(text, modality) swap_prefix(text, "response", modality)
to_package_names <- function(text, modality) swap_prefix(text, modality, "response")
rename_columns <- function(frame, map) {
    hit <- match(names(frame), names(map))
    names(frame)[!is.na(hit)] <- unname(map[hit[!is.na(hit)]])
    frame
}

SIDES <- c(left = -1L, both = 0L, right = 1L)
# The user-facing name "spa" is the default saddlepoint route of each mechanism.
SPA_DEFAULT <- c(crt = "crt_spa_fast", permutations = "rpt_spa_fast")

formula_text <- function(x) {
    if (inherits(x, "formula")) return(paste(deparse(x, width.cutoff = 500L), collapse = ""))
    if (is.list(x)) return(lapply(x, formula_text))
    x
}
text_formula <- function(x) {
    if (is.list(x)) {
        is_formula <- names(x) %in% "formula_object"
        x[is_formula] <- lapply(x[is_formula], text_formula)
        return(x)
    }
    stats::as.formula(x, env = globalenv())
}

rows_sorted <- function(table) {
    table <- as.data.frame(table, stringsAsFactors = FALSE)
    if (nrow(table) == 0L) return(table)
    table <- table[do.call(order, c(unname(as.list(table)), list(method = "radix"))), , drop = FALSE]
    rownames(table) <- NULL
    table
}
same <- function(a, b) isTRUE(all.equal(a, b, check.attributes = FALSE))
must_follow <- function(derived, held, what) {
    if (!same(derived, held)) fmt_stop("The object's `", what, "` does not follow from what the files hold.")
}

package_provenance <- function() {
    list(package = "sceptre", version = as.character(utils::packageVersion("sceptre")), language = "R")
}

# The data of a sceptre_object as the SingleCellExperiment fmt_write_screen takes.
object_screen <- function(object, response_names = NULL, id = NULL, response_modality = DEFAULT_RESPONSE_MODALITY) {
    ns <- sceptre_ns()
    screen <- ns$file_screen(object)
    meta <- S4Vectors::metadata(screen)
    grna_table <- as.data.frame(meta$grna_target_data_frame, stringsAsFactors = FALSE)
    # a gRNA under two targets has no place in the gRNA table (S8)
    if (anyDuplicated(grna_table$grna_id) > 0L) fmt_stop("A gRNA listed under two targets cannot be written in this format.")
    grna <- SingleCellExperiment::altExp(screen, "gRNA")
    grna_table <- grna_table[match(rownames(grna), grna_table$grna_id), , drop = FALSE]
    grna_table$grna_group <- NULL
    held <- as.list(SummarizedExperiment::colData(grna))
    if (length(held) > 0L) {
        must_follow(top_grna(SummarizedExperiment::assay(grna, "counts")), held[c("max_grna", "max_grna_frac_umis")], "top gRNA per cell")
    }
    SummarizedExperiment::colData(grna) <- NULL
    SummarizedExperiment::rowData(grna) <- S4Vectors::DataFrame(
        rename_columns(grna_table[setdiff(names(grna_table), "grna_id")], GRNA_COORDINATES),
        row.names = grna_table$grna_id, check.names = FALSE
    )
    SingleCellExperiment::altExps(screen) <- list()
    if (!is.null(response_names)) SummarizedExperiment::rowData(screen)$response_name <- as.character(response_names)
    SingleCellExperiment::altExp(screen, "gRNA") <- grna
    names(SummarizedExperiment::colData(screen)) <- to_file_names(names(SummarizedExperiment::colData(screen)), response_modality)
    S4Vectors::metadata(screen) <- list(sceptre = list(
        moi = if (isTRUE(meta$low_moi)) "low" else "high", id = id, response_modality = response_modality
    ))
    screen
}

write_sceptre_screen <- function(object, file, response_names = NULL, id = attr(object, "screen_id"),
                                 response_modality = DEFAULT_RESPONSE_MODALITY) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    fmt_write_screen(object_screen(object, response_names, id, response_modality), file)
}

called <- function(object) {
    done <- field(object, "functs_called")
    names(STAGE_FUNCTIONS)[STAGE_FUNCTIONS %in% names(done)[done]]
}

unit_columns <- function(grna_group, grna_table, strategy) {
    if (identical(strategy, "union")) return(data.frame(grna_target = as.character(grna_group)))
    data.frame(
        grna_target = grna_table$grna_target[match(grna_group, grna_table$grna_id)],
        grna_id = as.character(grna_group)
    )
}
unit_group <- function(table, strategy) if (identical(strategy, "union")) table$grna_target else table$grna_id

# Write the analysis of a sceptre_object into a new folder, as the stages
# would have written it one after another.
write_sceptre_analysis <- function(object, screen_file, directory, provenance = attr(object, "format_provenance"),
                                   response_modality = NULL) {
    ns <- sceptre_ns()
    done <- called(object)
    sce <- fmt_read_screen(screen_file, response_modality)
    modality <- S4Vectors::metadata(sce)$sceptre$response_modality
    cells <- colnames(sce)
    response_ids <- rownames(sce)
    grna_ids <- rownames(SingleCellExperiment::altExp(sce, "gRNA"))
    must_follow(colnames(ns$file_screen(object)), cells, "cells")
    must_follow(rownames(ns$file_screen(object)), response_ids, "responses")
    fmt_init_analysis(directory, screen_link(sce, screen_file, directory), modality)
    stage_provenance <- function(stage) if (is.null(provenance[[stage]])) package_provenance() else provenance[[stage]]
    write <- function(stage, record, tables = list()) {
        fmt_write_stage(directory, stage, record, tables[!vapply(tables, is.null, logical(1))], stage_provenance(stage))
    }
    low_moi <- isTRUE(field(object, "low_moi"))

    if ("parameters" %in% done) {
        strategy <- field(object, "grna_integration_strategy")
        grna_table <- as.data.frame(field(object, "grna_target_data_frame"), stringsAsFactors = FALSE)
        must_follow(grna_groups(grna_table, strategy), grna_table$grna_group, "grna_group")
        mechanism <- if (field(object, "run_permutations")) "permutations" else "crt"
        approximation <- field(object, "resampling_approximation")
        if (identical(approximation, SPA_DEFAULT[[mechanism]])) approximation <- "spa"
        # in the file, a covariate carries its modality's name
        formula <- to_file_names(formula_text(field(object, "formula_object")), modality)
        if (!formula_in_subset(formula)) {
            fmt_stop("The formula ", formula, " cannot be written: the format takes additive terms only, each a covariate name, log(x) or log(x + 1).")
        }
        design <- field(object, "covariate_matrix")
        rebuilt <- ns$convert_covariate_df_to_design_matrix(field(object, "covariate_data_frame"), field(object, "formula_object"))
        if (!identical(unname(design), unname(rebuilt))) fmt_stop("The object's design matrix is not the one its formula and covariates give.")
        pairs <- lapply(c(discovery = "discovery_pairs", positive_control = "positive_control_pairs"), function(name) {
            held <- as.data.frame(field(object, name), stringsAsFactors = FALSE)
            table <- unique(data.frame(grna_target = as.character(held$grna_target), response_id = as.character(held$response_id)))
            must_follow(rows_sorted(expand_pairs(table, grna_table)), rows_sorted(held[c("grna_target", "response_id", "grna_group")]), name)
            table
        })
        write("parameters", list(
            settings = list(
                formula = formula,
                side = names(SIDES)[match(field(object, "side_code"), SIDES)],
                resampling_mechanism = mechanism,
                resampling_approximation = approximation,
                control_group = if (field(object, "control_group_complement")) "complement" else "nt_cells",
                grna_integration_strategy = strategy,
                B1 = field(object, "B1"), B2 = field(object, "B2"), B3 = field(object, "B3"),
                multiple_testing_method = field(object, "multiple_testing_method"),
                multiple_testing_alpha = field(object, "multiple_testing_alpha"),
                seed = field(object, "seed"),
                response_fit_method = field(object, "response_fit_method"),
                grna_fit_method = field(object, "grna_fit_method")
            ),
            design = list(terms = I(to_file_names(colnames(design), modality)), sha256 = design_sha256(design))
        ), list("pairs/discovery" = pairs$discovery, "pairs/positive_control" = pairs$positive_control))
    }

    if ("assignment" %in% done) {
        initial <- field(object, "initial_grna_assignment_list")
        must_follow(names(initial), grna_ids, "the order of initial_grna_assignment_list")
        must_follow(lapply(initial, sort), initial, "the order of cells in initial_grna_assignment_list")
        must_follow(group_cells(initial, grna_table), field(object, "grna_assignments_raw"), "grna_assignments_raw")
        must_follow(mean(lengths(initial)), field(object, "mean_cells_per_grna"), "mean_cells_per_grna")
        held <- field(object, "grnas_per_cell")
        if (length(held) > 0L) must_follow(grnas_per_cell(initial, length(cells)), held, "grnas_per_cell")
        # cells the method itself flags as not carrying exactly one gRNA, beyond
        # what the assignment table shows (the "maximum" method does this)
        flagged <- setdiff(field(object, "cells_w_zero_or_twoplus_grnas"), if (low_moi) not_one_grna(initial, length(cells)))
        if (!low_moi && length(flagged) > 0L) fmt_stop("Cells are flagged in a high-MOI analysis.")
        hyperparameters <- formula_text(field(object, "grna_assignment_hyperparameters"))
        write("assignment", list(settings = list(
            method = field(object, "grna_assignment_method"),
            hyperparameters = if (length(hyperparameters) == 0L) stats::setNames(list(), character()) else hyperparameters
        )), list(
            assignments = data.frame(
                grna_id = rep(names(initial), lengths(initial)), cell_id = cells[unlist(initial, use.names = FALSE)]
            ),
            assignment_flagged_cells = if (length(flagged) > 0L) data.frame(cell_id = cells[flagged])
        ))
    }

    if ("qc" %in% done) {
        in_use <- field(object, "cells_in_use")
        nt_cells_control <- !field(object, "control_group_complement")
        must_follow(
            qc_assignments(field(object, "grna_assignments_raw"), in_use, nt_cells_control),
            field(object, "grna_assignments"), "grna_assignments"
        )
        pair_qc <- lapply(c(discovery = "discovery", positive_control = "positive_control"), function(type) {
            info <- as.data.frame(field(object, paste0(type, "_pairs_with_info")), stringsAsFactors = FALSE)
            table <- cbind(unit_columns(info$grna_group, grna_table, strategy), response_id = as.character(info$response_id), pass_qc = info$pass_qc)
            if (nrow(table) > 0L) must_follow(n_pairs_passing(table, strategy), field(object, paste0("n_ok_", type, "_pairs")), paste0("n_ok_", type, "_pairs"))
            table
        })
        removed <- field(object, "cell_removal_metrics")
        write("qc", list(
            settings = c(field(object, "cellwise_qc_thresholds"), list(pairwise_qc = field(object, "pairwise_qc"), B3 = field(object, "B3"))),
            summaries = list(cells_removed = as.list(stats::setNames(removed, sub("^n_cells_rm_", "", names(removed)))))
        ), list(
            cells = data.frame(cell_id = cells, in_use = seq_along(cells) %in% in_use),
            "pair_qc/discovery" = pair_qc$discovery, "pair_qc/positive_control" = pair_qc$positive_control
        ))
    }

    fits <- ns$response_fits(object)
    fits_stage <- intersect(c("calibration", "power", "discovery"), done)[1]
    fit_record <- function(stage) {
        if (length(fits) == 0L || !identical(stage, fits_stage)) return(NULL)
        list(file = paste0("fits/", stage, ".parquet"), identity = fit_identity(object, sce))
    }
    fit_table <- function(stage) {
        if (length(fits) == 0L || !identical(stage, fits_stage)) return(NULL)
        terms <- colnames(field(object, "covariate_matrix"))
        table <- data.frame(response_id = names(fits), theta = vapply(fits, function(f) f$theta, numeric(1)))
        table$coefficients <- lapply(fits, function(f) {
            must_follow(names(f$fitted_coefs), terms, "the terms of a saved fit")
            unname(f$fitted_coefs)
        })
        table[c("response_id", "coefficients", "theta")]
    }
    result_table <- function(stage) {
        result <- as.data.frame(ns$analysis_result(object, stage), stringsAsFactors = FALSE)
        if (!is.null(result$significant)) {
            must_follow(add_significant(result, significance_threshold(result))$significant, result$significant, "significant")
            result$significant <- NULL
        }
        # stored only where it does not follow from `fold_change` (T5)
        if (log_2_follows(result)) result$log_2_fold_change <- NULL
        result[setdiff(names(result), resample_columns(result))]
    }
    resamples <- function(stage) resamples_table(as.data.frame(ns$analysis_result(object, stage), stringsAsFactors = FALSE))

    if ("calibration" %in% done) {
        result <- as.data.frame(ns$analysis_result(object, "calibration"), stringsAsFactors = FALSE)
        negative <- as.data.frame(field(object, "negative_control_pairs"), stringsAsFactors = FALSE)
        # a control group is named in the column that names the tested unit:
        # `grna_target` under union, `grna_id` otherwise
        must_follow(
            rows_sorted(data.frame(response_id = result$response_id, grna_group = unit_group(result, strategy), pass_qc = result$pass_qc)),
            rows_sorted(data.frame(response_id = negative$response_id, grna_group = negative$grna_group, pass_qc = negative$pass_qc)),
            "negative_control_pairs"
        )
        # today a control group's members exist only in its name (T4)
        groups <- unique(as.character(negative$grna_group))
        members <- strsplit(groups, "&", fixed = TRUE)
        write("calibration", list(
            settings = list(
                calibration_group_size = field(object, "calibration_group_size"),
                n_calibration_pairs = field(object, "n_calibration_pairs"), B3 = field(object, "B3")
            ),
            significance_threshold = significance_threshold(result),
            fits = fit_record("calibration")
        ), list(
            negative_control_groups = data.frame(grna_group = rep(groups, lengths(members)), grna_id = unlist(members)),
            "results/calibration" = result_table("calibration"), "results/calibration_resamples" = resamples("calibration"),
            "fits/calibration" = fit_table("calibration")
        ))
    }
    if ("power" %in% done) {
        write("power", list(settings = stats::setNames(list(), character()), fits = fit_record("power")),
            list("results/power" = result_table("power"), "results/power_resamples" = resamples("power"), "fits/power" = fit_table("power")))
    }
    if ("discovery" %in% done) {
        write("discovery", list(
            settings = stats::setNames(list(), character()),
            significance_threshold = significance_threshold(as.data.frame(ns$analysis_result(object, "discovery"))),
            fits = fit_record("discovery")
        ), list(
            "results/discovery" = result_table("discovery"), "results/discovery_resamples" = resamples("discovery"),
            "fits/discovery" = fit_table("discovery")
        ))
    }
    invisible(directory)
}

# What makes a saved fit valid (T6).
fit_identity <- function(object, sce) {
    list(
        screen_id = S4Vectors::metadata(sce)$sceptre$id,
        formula = to_file_names(formula_text(field(object, "formula_object")), S4Vectors::metadata(sce)$sceptre$response_modality),
        cells_in_use_sha256 = sha256_ids(colnames(sce)[field(object, "cells_in_use")]),
        response_fit_method = field(object, "response_fit_method"),
        control_group = if (field(object, "control_group_complement")) "complement" else "nt_cells"
    )
}

# Read a data file, and optionally one analysis of it, as a sceptre_object.
# Without `directory` the object is as after import_data(). The response
# modality is the one the analysis names; without an analysis, `response_modality`
# or the file's only modality besides `grna`.
read_sceptre_object <- function(screen_file, directory = NULL, response_modality = NULL) {
    ns <- sceptre_ns()
    if (!is.null(directory)) response_modality <- fmt_response_modality(directory)
    sce <- fmt_read_screen(screen_file, response_modality)
    modality <- S4Vectors::metadata(sce)$sceptre$response_modality
    grna <- SingleCellExperiment::altExp(sce, "gRNA")
    grna_table <- data.frame(
        grna_id = rownames(grna),
        rename_columns(
            as.data.frame(SummarizedExperiment::rowData(grna), optional = TRUE),
            stats::setNames(names(GRNA_COORDINATES), GRNA_COORDINATES)
        ),
        stringsAsFactors = FALSE
    )
    rownames(grna_table) <- NULL
    cells <- colnames(sce)
    response_ids <- rownames(sce)
    low_moi <- identical(S4Vectors::metadata(sce)$sceptre$moi, "low")
    covariates <- SummarizedExperiment::colData(sce)
    names(covariates) <- to_package_names(names(covariates), modality)
    file_screen <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = ns$as_row_counts(SummarizedExperiment::assay(sce, "counts"))),
        colData = covariates,
        altExps = list(gRNA = SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = ns$as_row_counts(SummarizedExperiment::assay(grna, "counts"))),
            colData = S4Vectors::DataFrame(top_grna(SummarizedExperiment::assay(grna, "counts")), row.names = cells)
        )),
        metadata = list(low_moi = low_moi, grna_target_data_frame = grna_table, covariate_row_names = cells)
    )
    object <- ns$object_from_file_screen(file_screen)
    # what import_data() sets beyond the data
    field(object, "covariate_names") <- sort(names(field(object, "covariate_data_frame")))
    field(object, "nuclear") <- FALSE
    field(object, "nf_pipeline") <- FALSE
    field(object, "elements_to_analyze") <- NA_character_
    done <- character()
    provenance <- list()
    if (!is.null(directory)) {
        analysis <- fmt_read_analysis(directory)
        check_link(analysis$screen, sce, directory)
        done <- names(analysis$stages)
        tables <- analysis$tables
        records <- analysis$stages
        provenance <- lapply(records, function(r) r$provenance)
    }
    set <- function(name, value) field(object, name) <<- value

    if ("parameters" %in% done) {
        settings <- records$parameters$settings
        strategy <- settings$grna_integration_strategy
        grna_table$grna_group <- grna_groups(grna_table, strategy)
        approximation <- settings$resampling_approximation
        if (identical(approximation, "spa")) approximation <- SPA_DEFAULT[[settings$resampling_mechanism]]
        set("grna_target_data_frame", grna_table)
        set("formula_object", text_formula(to_package_names(settings$formula, modality)))
        set("side_code", SIDES[[settings$side]])
        set("run_permutations", identical(settings$resampling_mechanism, "permutations"))
        set("resampling_approximation", approximation)
        set("control_group_complement", identical(settings$control_group, "complement"))
        set("grna_integration_strategy", strategy)
        for (name in c("B1", "B2", "B3", "seed")) set(name, as.integer(settings[[name]]))
        set("multiple_testing_method", settings$multiple_testing_method)
        set("multiple_testing_alpha", as.double(settings$multiple_testing_alpha))
        set("response_fit_method", settings$response_fit_method)
        set("grna_fit_method", settings$grna_fit_method)
        design <- ns$convert_covariate_df_to_design_matrix(field(object, "covariate_data_frame"), field(object, "formula_object"))
        if (!identical(to_file_names(colnames(design), modality), as.character(records$parameters$design$terms)) ||
            !identical(design_sha256(design), records$parameters$design$sha256)) {
            fmt_stop("The design matrix rebuilt from the formula is not the one the analysis was run with.")
        }
        set("covariate_matrix", design)
        for (type in c("discovery", "positive_control")) {
            pairs <- expand_pairs(tables[[paste0("pairs/", type)]], grna_table)
            rownames(pairs) <- NULL
            set(paste0(type, "_pairs"), pairs)
            # counted as the user gave them, before any expansion by gRNA
            set(paste0("n_", type, "_pairs"), nrow(tables[[paste0("pairs/", type)]]))
        }
    }

    if ("assignment" %in% done) {
        settings <- records$assignment$settings
        set("grna_assignment_method", settings$method)
        set("grna_assignment_hyperparameters", text_formula(as.list(settings$hyperparameters)))
        initial <- cells_by_grna(tables$assignments, rownames(grna), cells)
        raw <- group_cells(initial, grna_table)
        set("initial_grna_assignment_list", initial)
        set("grna_assignments_raw", raw)
        set("mean_cells_per_grna", mean(lengths(initial)))
        # the "maximum" method does not count gRNAs per cell
        if (!identical(settings$method, "maximum")) set("grnas_per_cell", grnas_per_cell(initial, length(cells)))
        if (low_moi) {
            flagged <- match(tables$assignment_flagged_cells$cell_id, cells)
            set("cells_w_zero_or_twoplus_grnas", sort(union(not_one_grna(initial, length(cells)), flagged)))
        }
    }

    if ("qc" %in% done) {
        settings <- records$qc$settings
        set("cellwise_qc_thresholds", settings[c("response_n_umis_range", "response_n_nonzero_range", "p_mito_threshold")])
        set("pairwise_qc", isTRUE(settings$pairwise_qc))
        set("B3", as.integer(settings$B3))
        removed <- unlist(records$qc$summaries$cells_removed)
        set("cell_removal_metrics", stats::setNames(as.integer(removed), paste0("n_cells_rm_", names(removed))))
        table_cells <- tables$cells[match(cells, tables$cells$cell_id), , drop = FALSE]
        in_use <- which(table_cells$in_use)
        set("cells_in_use", in_use)
        set("grna_assignments", qc_assignments(raw, in_use, identical(records$parameters$settings$control_group, "nt_cells")))
        for (type in c("discovery", "positive_control")) {
            pair_qc <- tables[[paste0("pair_qc/", type)]]
            info <- data.frame(response_id = pair_qc$response_id, grna_group = unit_group(pair_qc, strategy), pass_qc = pair_qc$pass_qc)
            # the order a package tests in: by response, in the data file's order
            info <- info[order(match(info$response_id, response_ids)), , drop = FALSE]
            rownames(info) <- NULL
            set(paste0(type, "_pairs_with_info"), data.table::as.data.table(info))
            if (nrow(info) > 0L) set(paste0("n_ok_", type, "_pairs"), as.integer(n_pairs_passing(pair_qc, strategy)))
        }
        # per response, the cells of the control pool that express it: a pass
        # over the counts, done with the package's own routine
        object <- restore_pool_nonzero(object)
    }

    result <- function(stage) {
        table <- add_log_2_fold_change(tables[[paste0("results/", stage)]])
        table <- add_resamples(table, tables[[paste0("results/", stage, "_resamples")]])
        if (stage != "power") table <- add_significant(table, records[[stage]]$significance_threshold)
        # as sceptre returns results: by p-value, then response
        table <- table[do.call(order, c(list(table$p_value, table$response_id), unname(as.list(table[unit_key(table)])), list(method = "radix"))), , drop = FALSE]
        rownames(table) <- NULL
        data.table::as.data.table(table)
    }
    set_result <- function(stage) {
        setter <- get("analysis_result<-", envir = ns)
        object <<- setter(object, stage, value = result(stage))
    }
    if ("calibration" %in% done) {
        settings <- records$calibration$settings
        set("calibration_group_size", as.integer(settings$calibration_group_size))
        set("n_calibration_pairs", as.integer(settings$n_calibration_pairs))
        set("B3", as.integer(settings$B3))
        calibration <- tables[["results/calibration"]]
        set("negative_control_pairs", data.table::as.data.table(data.frame(
            response_id = calibration$response_id, grna_group = unit_group(calibration, strategy), pass_qc = calibration$pass_qc
        )))
        set_result("calibration")
    }
    if ("power" %in% done) set_result("power")
    if ("discovery" %in% done) set_result("discovery")

    # a saved fit is reused only if it was made as this analysis would make it (T6)
    for (stage in intersect(c("calibration", "power", "discovery"), done)) {
        saved <- records[[stage]]$fits
        if (is.null(saved)) next
        reason <- fit_mismatch(saved$identity, fit_identity(object, sce), records[[stage]]$provenance$language)
        if (!is.null(reason)) {
            message("Saved fits in ", saved$file, " are not reused: ", reason, ".")
            next
        }
        table <- tables[[paste0("fits/", stage)]]
        terms <- colnames(field(object, "covariate_matrix"))
        fits <- lapply(stats::setNames(seq_len(nrow(table)), table$response_id), function(i) {
            list(fitted_coefs = stats::setNames(table$coefficients[[i]], terms), theta = table$theta[i])
        })
        setter <- get("response_fits<-", envir = ns)
        object <- setter(object, value = c(ns$response_fits(object), fits))
    }

    status <- stats::setNames(c(TRUE, STAGE_FUNCTIONS %in% STAGE_FUNCTIONS[done]), c("import_data", STAGE_FUNCTIONS))
    field(object, "functs_called") <- status
    field(object, "last_function_called") <- if (length(done) == 0L) "import_data" else STAGE_FUNCTIONS[[utils::tail(done, 1L)]]
    attr(object, "format_provenance") <- provenance
    attr(object, "screen_id") <- S4Vectors::metadata(sce)$sceptre$id
    object
}

fit_mismatch <- function(saved, current, language) {
    for (name in names(current)) {
        if (!identical(saved[[name]], current[[name]])) return(paste0("their `", name, "` differs"))
    }
    method <- current$response_fit_method
    if (!"R" %in% FIT_METHODS[[method]] || !language %in% FIT_METHODS[[method]]) {
        return(paste0("the method '", method, "' is not shared between ", language, " and R"))
    }
    NULL
}

# Per response, how many cells of the control pool express it. sceptre computes
# this during pairwise QC; it is recomputed here with sceptre's own routine.
restore_pool_nonzero <- function(object) {
    ns <- sceptre_ns()
    info <- rbind(
        as.data.frame(field(object, "discovery_pairs_with_info")),
        as.data.frame(field(object, "positive_control_pairs_with_info"))
    )
    if (nrow(info) == 0L) return(object)
    response_ids <- rownames(ns$get_response_matrix(object))
    groups <- names(field(object, "grna_assignments")$grna_group_idxs)
    # as sceptre calls it: pairs ordered by response
    info <- info[order(match(info$response_id, response_ids)), , drop = FALSE]
    margins <- ns$pair_margins(object, match(info$response_id, response_ids), match(info$grna_group, groups))
    field(object, "n_nonzero_tot_vector") <- margins$pool_nonzero
    object
}

# Compare two sceptre_objects field by field. Tables are compared with their
# rows sorted and the fit cache by response, since the files keep neither order.
compare_objects <- function(a, b) {
    ns <- sceptre_ns()
    fields <- setdiff(ns$object_field_names, "response_precomputations")
    plain <- function(x, name) {
        if (inherits(x, "formula")) return(paste(deparse(x, width.cutoff = 500L), collapse = ""))
        if (is.data.frame(x)) return(rows_sorted(x))
        if (methods::is(x, "sparseMatrix")) return(methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix"))
        if (name == "cells_w_zero_or_twoplus_grnas") return(sort(x))
        x
    }
    verdict <- vapply(fields, function(name) {
        x <- plain(field(a, name), name)
        y <- plain(field(b, name), name)
        if (identical(x, y)) return("identical")
        difference <- all.equal(x, y, check.attributes = FALSE)
        if (isTRUE(difference)) "equal" else paste(difference, collapse = "; ")
    }, character(1))
    by_name <- function(fits) fits[sort(names(fits))]
    fits <- if (identical(by_name(ns$response_fits(a)), by_name(ns$response_fits(b)))) {
        "identical"
    } else {
        difference <- all.equal(by_name(ns$response_fits(a)), by_name(ns$response_fits(b)), check.attributes = FALSE)
        if (isTRUE(difference)) "equal" else paste(difference, collapse = "; ")
    }
    c(verdict, response_fits = fits)
}

# The files of a folder with their sizes.
file_listing <- function(directory) {
    files <- list.files(directory, recursive = TRUE)
    data.frame(file = files, kilobytes = round(file.size(file.path(directory, files)) / 1024, 1))
}

# A one-line summary of compare_objects().
comparison_summary <- function(verdict) {
    different <- names(verdict)[!verdict %in% c("identical", "equal")]
    paste0(
        sum(verdict == "identical"), " fields identical, ", sum(verdict == "equal"), " equal, ", length(different), " different",
        if (length(different) > 0L) paste0(": ", paste(different, collapse = ", ")) else ""
    )
}

# Whether two tables hold the same rows, in any order.
same_table <- function(a, b) identical(rows_sorted(a), rows_sorted(b))

# The integer width of each array of each count matrix in a data file.
count_types <- function(file) {
    h5 <- rhdf5::H5Fopen(file, flags = "H5F_ACC_RDONLY")
    on.exit(rhdf5::H5Fclose(h5), add = TRUE)
    modalities <- sort(unique(sub("^/mod/?", "", subset(rhdf5::h5ls(file, recursive = 2L), group == "/mod")$name)))
    paths <- as.vector(outer(c("data", "indices", "indptr"), modalities, function(name, modality) paste0("mod/", modality, "/X/", name)))
    vapply(stats::setNames(paths, sub("^mod/", "", sub("/X/", " ", paths))), function(path) {
        dataset <- rhdf5::H5Dopen(h5, path)
        on.exit(rhdf5::H5Dclose(dataset), add = TRUE)
        paste0(8L * rhdf5::H5Tget_size(rhdf5::H5Dget_type(dataset)), "-bit")
    }, character(1))
}

# ---------------------------------------------------------------- validation

# validate_screen() and validate_analysis() check files against the numbered
# rules of SPEC.md. They read the files directly, build no object, and report
# every problem they find: a data frame with one row per problem (`rule`,
# `where`, `message`), with no rows when the files conform. The Python library
# has the same two functions and reports the same rules.

# `readable` marks a problem that makes a file not conforming while a reader
# must still accept it (another whole-number type for the counts, S4.4).
problem_log <- function() {
    found <- list()
    record <- function(readable) function(rule, where, ...) {
        found[[length(found) + 1L]] <<- data.frame(rule = rule, where = where, message = paste0(...), readable = readable)
    }
    list(
        add = record(FALSE), add_readable = record(TRUE),
        all = function() {
            if (length(found) == 0L) return(data.frame(rule = character(), where = character(), message = character(), readable = logical()))
            out <- unique(do.call(rbind, found))
            out <- out[order(out$rule, out$where, out$message, method = "radix"), ]
            rownames(out) <- NULL
            out
        }
    )
}

MODALITY_NAME <- SCREEN_SPEC$modalities$name_pattern

# One word for what a set of problems means: "conforming" (none), "readable"
# (not conforming, and a reader must still accept the file) or "not readable".
conformance_level <- function(problems) {
    if (nrow(problems) == 0L) "conforming" else if (all(problems$readable)) "readable" else "not readable"
}

# The compression filters of HDF5 datasets.
h5_filters <- function(file, paths) {
    h5 <- rhdf5::H5Fopen(file, flags = "H5F_ACC_RDONLY")
    on.exit(rhdf5::H5Fclose(h5), add = TRUE)
    lapply(stats::setNames(paths, paths), function(path) {
        dataset <- rhdf5::H5Dopen(h5, path)
        on.exit(rhdf5::H5Dclose(dataset), add = TRUE)
        plist <- rhdf5::H5Dget_create_plist(dataset)
        on.exit(rhdf5::H5Pclose(plist), add = TRUE, after = FALSE)
        vapply(seq_len(rhdf5::H5Pget_nfilters(plist)), function(i) rhdf5::H5Pget_filter(plist, i)[[2]], character(1))
    })
}
# The storage type of a dataset, from its HDF5 type name.
STORAGE <- c(int32 = "^H5T_STD_I32", int64 = "^H5T_STD_I64", double = "^H5T_IEEE_F64")
is_stored_as <- function(dtype, kind) grepl(STORAGE[[kind]], dtype)
# The kind of a column of a cell or feature table (spec/screen.json, storage_kinds).
storage_kind <- function(encoding, dtype) {
    if (identical(encoding, "categorical")) return("categorical")
    if (identical(encoding, "string-array")) return("string")
    if (identical(encoding, "array") && is_stored_as(dtype, "double")) return("double")
    if (identical(encoding, "array") && is_stored_as(dtype, "int32")) return("int32")
    "other"
}

validate_screen <- function(file) {
    log <- problem_log()
    listing <- tryCatch(rhdf5::h5ls(file, all = TRUE), error = function(e) NULL)
    if (is.null(listing)) {
        log$add("S1.1", file, "is not an HDF5 file.")
        return(log$all())
    }
    attribute <- function(path, name) tryCatch(rhdf5::h5readAttributes(file, path)[[name]], error = function(e) NULL)
    has <- function(path) {
        path <- paste0("/", sub("^/", "", path))
        any(listing$group == dirname(path) & listing$name == basename(path)) || path == "/"
    }
    children <- function(path) listing[listing$group == paste0("/", sub("^/", "", path)), , drop = FALSE]
    read <- function(path) as.vector(rhdf5::h5read(file, path, bit64conversion = "double"))
    if (!identical(as.vector(attribute("/", "encoding-type")), SCREEN_SPEC$root[["encoding-type"]])) log$add("S1.1", file, "The file is not marked as MuData.")
    if (!has("mod")) {
        log$add("S2.1", file, "The file has no modalities.")
        return(log$all())
    }
    in_file <- children("mod")$name
    if (!GRNA %in% in_file || length(in_file) < SCREEN_SPEC$modalities$minimum) log$add("S2.1", "mod", "A data file needs the modality 'grna' and at least one other.")
    for (name in in_file[!grepl(MODALITY_NAME, in_file)]) log$add("S2.2", paste0("mod/", name), "is not a valid modality name.")

    # file-level metadata
    where <- SCREEN_SPEC$metadata$path
    record <- children(where)
    entry <- function(name) record[record$name == name, , drop = FALSE]
    scalar_of <- function(name, class) nrow(entry(name)) == 1L && entry(name)$dclass == class && identical(trimws(entry(name)$dim), "( 0 )")
    entries <- SCREEN_SPEC$metadata$entries
    classes <- vapply(entries, function(e) if (e$kind == "integer") "INTEGER" else "STRING", character(1))
    if (!all(mapply(scalar_of, names(entries), classes))) {
        log$add("S10.1", where, "needs the scalar text entries `format`, `moi` and `id`.")
    } else if (is.null(format_parts(read(paste0(where, "/format"))))) {
        log$add("S10.1", paste0(where, "/format"), "is not a version of the form MAJOR.MINOR.")
    } else {
        if (!read(paste0(where, "/moi")) %in% unlist(entries$moi$values)) log$add("S10.2", paste0(where, "/moi"), "must be 'high' or 'low'.")
        if (!known_major(read(paste0(where, "/format")))) log$add("S10.3", paste0(where, "/format"), "is a format this validator does not know.")
        if (!nzchar(read(paste0(where, "/id")))) log$add("S10.1", paste0(where, "/id"), "is empty.")
    }

    index_of <- function(table) {
        name <- as.vector(attribute(table, "_index"))
        if (is.null(name) || !has(paste0(table, "/", name))) return(NULL)
        as.character(read(paste0(table, "/", name)))
    }
    # cells: the same, in the same order, in every modality and in the shared table
    cells <- if (GRNA %in% in_file) index_of("mod/grna/obs") else index_of(paste0("mod/", in_file[1], "/obs"))
    if (is.null(cells)) {
        log$add("S2.3", "mod", "A modality has no cell identifiers.")
        return(log$all())
    }
    if (anyDuplicated(cells) > 0L) log$add("S9.1", "obs", "Cell identifiers are not unique.")
    for (name in in_file) {
        if (!identical(index_of(paste0("mod/", name, "/obs")), cells)) log$add("S2.3", paste0("mod/", name, "/obs"), "does not hold the same cells in the same order as 'grna'.")
    }
    if (has("obs") && !identical(index_of("obs"), cells)) log$add("S2.3", "obs", "The shared cell table does not hold the same cells in the same order.")
    features <- lapply(stats::setNames(in_file, in_file), function(name) index_of(paste0("mod/", name, "/var")))
    if (anyDuplicated(unlist(features)) > 0L) log$add("S9.2", "mod", "Feature identifiers are not unique across the modalities.")
    if (any(grepl("\n", c(cells, unlist(features)), fixed = TRUE))) log$add("S9.4", "mod", "An identifier contains a line feed.")

    # the count matrix of every modality (S3, S4, S11)
    expect <- unlist(SCREEN_SPEC$counts$arrays)
    rules <- unlist(SCREEN_SPEC$counts$rules)
    for (name in in_file) {
        x <- paste0("mod/", name, "/X")
        n_features <- length(features[[name]])
        if (!identical(as.vector(attribute(x, "encoding-type")), SCREEN_SPEC$counts[["encoding-type"]]) ||
            !identical(as.numeric(attribute(x, "shape")), as.numeric(c(length(cells), n_features)))) {
            log$add("S3.1", x, "is not a column-compressed (csc_matrix) matrix of cells x features.")
            next
        }
        filters <- h5_filters(file, paste0(x, "/", names(expect)))
        stored <- children(x)
        for (array in names(expect)) {
            row <- stored[stored$name == array, , drop = FALSE]
            if (!is_stored_as(row$dtype, expect[[array]])) {
                # another whole-number type: not conforming, and a reader still accepts it (S4.4)
                add <- if (row$dclass %in% c("INTEGER", "FLOAT")) log$add_readable else log$add
                add(rules[[array]], paste0(x, "/", array), "is not stored as ", sub("int", "", expect[[array]]), "-bit signed integers.")
            }
            if (length(setdiff(filters[[paste0(x, "/", array)]], unlist(SCREEN_SPEC$counts$filters))) > 0L) log$add("S11.1", paste0(x, "/", array), "uses a filter other than gzip.")
        }
        data <- read(paste0(x, "/data"))
        indices <- read(paste0(x, "/indices"))
        indptr <- read(paste0(x, "/indptr"))
        if (length(indptr) != n_features + 1L || indptr[1] != 0 || is.unsorted(indptr) || utils::tail(indptr, 1) != length(data) || length(indices) != length(data)) {
            log$add("S3.2", x, "The column pointers do not describe the stored values.")
            next
        }
        if (length(indices) > 0L) {
            step <- diff(indices)
            starts <- indptr[indptr > 0 & indptr < length(indices)]   # where a new feature begins
            within <- if (length(starts) > 0L) step[-starts] else step
            if (any(indices < 0 | indices >= length(cells)) || any(within <= 0)) log$add("S3.3", x, "Within a feature, cell positions must be in range and strictly increasing.")
            if (any(data <= 0) || any(data != round(data))) log$add("S3.4", x, "Stored counts must be positive integers.")
        }
    }

    # column kinds (S7); returns each column's kind
    check_table <- function(table, cell_table) {
        allowed <- unlist(if (cell_table) SCREEN_SPEC$cell_tables$column_kinds else SCREEN_SPEC$feature_tables$column_kinds)
        index <- as.vector(attribute(table, "_index"))
        entries <- children(table)
        entries <- entries[entries$name != index, , drop = FALSE]
        kinds <- stats::setNames(character(nrow(entries)), entries$name)
        for (i in seq_len(nrow(entries))) {
            path <- paste0(table, "/", entries$name[i])
            kind <- storage_kind(as.vector(attribute(path, "encoding-type")), entries$dtype[i])
            kinds[i] <- kind
            if (!kind %in% allowed) {
                log$add(if (cell_table) "S7.1" else "S7.2", path, "is not one of the column kinds allowed here (", paste(allowed, collapse = ", "), ").")
            } else if (kind == "categorical") {
                categories <- read(paste0(path, "/categories"))
                codes <- read(paste0(path, "/codes"))
                if (!is.character(categories) || anyDuplicated(categories) > 0L || any(codes < -1 | codes >= length(categories))) {
                    log$add("S7.3", path, "is not a well-formed categorical.")
                }
            } else if (kind == "int32") {
                # rhdf5 returns the value -2^31 as NA, with a message
                if (anyNA(suppressMessages(read(path)))) log$add("S7.4", path, "An integer column holds a missing value; it must be stored as doubles.")
            }
            if (cell_table && !grepl(SCREEN_SPEC$cell_tables$column_name_pattern, entries$name[i])) {
                log$add("S5.4", path, "is not a valid covariate name (letters, digits, underscores and dots, starting with a letter).")
            }
        }
        kinds
    }
    if (has("obs")) {
        shared <- names(check_table("obs", TRUE))
        for (name in in_file) {
            for (column in shared[startsWith(shared, paste0(name, "_"))]) {
                log$add("S5.1", paste0("obs/", column), "A shared column must not start with a modality's name and an underscore.")
            }
        }
    }
    for (name in in_file) {
        table <- paste0("mod/", name, "/obs")
        kinds <- check_table(table, TRUE)
        # the columns every modality's cell table has (S6)
        for (column in SCREEN_SPEC$cell_tables$every_modality) {
            if (!column$name %in% names(kinds)) {
                if (isTRUE(column$required)) log$add("S6.1", table, "has no column `", column$name, "`.")
            } else if (!kinds[[column$name]] %in% unlist(column$kinds)) {
                log$add("S6.2", paste0(table, "/", column$name), "must be a number.")
            }
        }
        check_table(paste0("mod/", name, "/var"), FALSE)
    }

    # the gRNA table (S8)
    if (GRNA %in% in_file) {
        table <- "mod/grna/var"
        columns <- children(table)
        text <- function(name) {
            path <- paste0(table, "/", name)
            encoding <- as.vector(attribute(path, "encoding-type"))
            if (identical(encoding, "string-array")) return(as.character(read(path)))
            codes <- read(paste0(path, "/codes"))
            ifelse(codes < 0, NA_character_, as.character(read(paste0(path, "/categories")))[codes + 1L])
        }
        for (column in SCREEN_SPEC$grna_table_forbidden) {
            if (column$name %in% columns$name) log$add(column$rule, paste0(table, "/", column$name), "must not be stored.")
        }
        for (column in SCREEN_SPEC$grna_table) {
            row <- columns[columns$name == column$name, , drop = FALSE]
            path <- paste0(table, "/", column$name)
            if (nrow(row) == 0L) {
                if (isTRUE(column$required)) log$add(column$rule, table, "The gRNA table has no `", column$name, "`.")
                next
            }
            kind <- storage_kind(as.vector(attribute(path, "encoding-type")), row$dtype)
            if (!kind %in% unlist(column$kinds)) {
                log$add(column$rule, path, "must be ", if ("string" %in% unlist(column$kinds)) "text" else "a number", if (isTRUE(column$required)) " with no missing values." else ".")
            } else if (isTRUE(column$required)) {
                values <- text(column$name)
                if (anyNA(values) || any(!nzchar(values))) log$add(column$rule, path, "must be text with no missing values.")
            }
        }
    }
    log$all()
}

# What validate_analysis() needs from a data file, read without building an object.
screen_identity <- function(file, response_modality) {
    read <- function(path) as.character(as.vector(rhdf5::h5read(file, path)))
    index <- function(table) read(paste0(table, "/", as.vector(rhdf5::h5readAttributes(file, table)[["_index"]])))
    modalities <- subset(rhdf5::h5ls(file, recursive = 2L), group == "/mod")$name
    out <- list(modalities = modalities, id = read("uns/sceptre/id"), cells = index("mod/grna/obs"), grnas = index("mod/grna/var"))
    target <- "mod/grna/var/grna_target"
    out$targets <- if (identical(as.vector(rhdf5::h5readAttributes(file, target)[["encoding-type"]]), "categorical")) {
        read(paste0(target, "/categories"))[as.integer(rhdf5::h5read(file, paste0(target, "/codes"))) + 1L]
    } else {
        read(target)
    }
    if (response_modality %in% modalities) out$responses <- index(paste0("mod/", response_modality, "/var"))
    out
}

json_validator <- local({
    cache <- list()
    function(name) {
        if (is.null(cache[[name]])) cache[[name]] <<- jsonvalidate::json_validator(file.path(SPEC_DIR, name), engine = "ajv")
        cache[[name]]
    }
})
# Problems a JSON file has against a schema, as text; character() when it conforms.
schema_problems <- function(path, schema) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    parsed <- tryCatch(jsonlite::fromJSON(text, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(parsed)) return("is not valid JSON.")
    result <- json_validator(schema)(text, verbose = TRUE, greedy = TRUE)
    if (isTRUE(result)) return(character())
    errors <- attr(result, "errors")
    unique(paste0(ifelse(nzchar(errors$instancePath), errors$instancePath, "/"), " ", errors$message))
}

# The tables a stage record lists, read without trusting the rest of the record.
listed_paths <- function(path) {
    files <- tryCatch(read_json(path)$files, error = function(e) NULL)
    paths <- unlist(lapply(files, function(f) if (is.list(f)) f$path))
    sub("\\.parquet$", "", as.character(paths))
}

column_kind <- function(column) {
    if (is.character(column) || is.factor(column)) return("string")
    if (is.logical(column)) return("boolean")
    if (is.double(column)) return("double")
    if (is.list(column) && all(vapply(column, function(v) is.null(v) || is.double(v), logical(1)))) return("list<double>")
    "other"
}

validate_analysis <- function(directory, screen_file = NULL) {
    log <- problem_log()
    path <- file.path(directory, "analysis.json")
    if (!file.exists(path)) {
        log$add("A1.1", "analysis.json", "is missing.")
        return(log$all())
    }
    faults <- schema_problems(path, "analysis.schema.json")
    for (fault in faults) log$add("A1.1", "analysis.json", fault)
    analysis <- if (length(faults) == 0L) read_json(path) else NULL
    if (!is.null(analysis) && !known_major(analysis$format)) log$add("A1.1", "analysis.json", "is in a format this validator does not know.")

    # stage records
    present <- sub("\\.json$", "", list.files(file.path(directory, "stages"), pattern = "\\.json$"))
    for (name in setdiff(present, STAGES)) log$add("A2.1", paste0("stages/", name, ".json"), "is not named for a stage.")
    records <- list()
    unchecked <- character()
    for (stage in intersect(STAGES, present)) {
        where <- paste0("stages/", stage, ".json")
        faults <- schema_problems(stage_file(directory, stage), "stage.schema.json")
        for (fault in faults) log$add("A2.1", where, fault)
        if (length(faults) > 0L) {
            # its files are still accounted for, so that they are not also reported as stale
            unchecked <- c(unchecked, listed_paths(stage_file(directory, stage)))
            next
        }
        record <- read_json(stage_file(directory, stage))
        if (!identical(record$stage, stage)) {
            log$add("A2.1", where, "says it is the record of the stage '", record$stage, "'.")
            unchecked <- c(unchecked, listed_paths(stage_file(directory, stage)))
            next
        }
        records[[stage]] <- record
    }
    for (stage in intersect(STAGES, present)) {
        missing <- setdiff(names(STAGE_RANK)[STAGE_RANK < STAGE_RANK[[stage]]], present)
        if (length(missing) > 0L) log$add("A2.2", paste0("stages/", stage, ".json"), "The stage needs the stage(s) ", paste(missing, collapse = ", "), " first.")
    }

    # the files each record lists
    listed <- character()
    tables <- list()
    for (stage in names(records)) {
        record <- records[[stage]]
        where <- paste0("stages/", stage, ".json")
        paths <- vapply(record$files, function(f) f$path, character(1))
        names_listed <- sub("\\.parquet$", "", paths)
        listed <- c(listed, names_listed)
        for (file in record$files) {
            name <- sub("\\.parquet$", "", file$path)
            known <- TABLES[[name]]
            if (!is.null(known) && known$stage != stage) log$add("A5.3", where, "lists ", file$path, ", which the stage '", known$stage, "' writes.")
            full <- file.path(directory, file$path)
            if (!file.exists(full)) {
                log$add("A5.1", file$path, "is listed by the stage '", stage, "' and is missing.")
                next
            }
            rows <- tryCatch(table_rows(full), error = function(e) NA)
            if (is.na(rows)) {
                log$add("A5.1", file$path, "cannot be read as Parquet.")
                next
            }
            if (rows != file$rows) log$add("A5.2", file$path, "has ", rows, " rows; the stage '", stage, "' wrote ", file$rows, ".")
            if (!is.null(known)) tables[[name]] <- full
        }
        required <- names(TABLES)[vapply(TABLES, function(t) t$stage == stage && isTRUE(t$required), logical(1))]
        for (name in setdiff(required, names_listed)) log$add("A5.4", where, "does not list the table ", name, ".")
        # saved fits and their record go together (T6)
        fit_table <- paste0("fits/", stage)
        if (!is.null(record$fits) && !identical(sub("\\.parquet$", "", record$fits$file), fit_table)) log$add("T6.1", where, "`fits` names a file other than ", fit_table, ".parquet.")
        if (is.null(record$fits) != !(fit_table %in% names_listed)) log$add("T6.1", where, "Saved fits and the record of how they were made must come together.")
    }
    for (name in setdiff(names(TABLES), c(listed, unchecked))) {
        if (file.exists(file.path(directory, paste0(name, ".parquet")))) log$add("A4.1", paste0(name, ".parquet"), "is in the folder and no stage record lists it.")
    }

    # each table against its definition
    strategy <- records$parameters$settings$grna_integration_strategy
    contents <- list()
    for (name in names(tables)) {
        definition <- TABLES[[name]]$definition
        where <- paste0(name, ".parquet")
        parts <- lapply(table_files(tables[[name]]), function(file) {
            part <- as.data.frame(nanoparquet::read_parquet(file), stringsAsFactors = FALSE)
            class(part) <- "data.frame"
            groups <- nanoparquet::read_parquet_metadata(file)$row_groups$num_rows
            if (any(groups > ROW_GROUP_ROWS)) log$add("T0.3", where, "has a row group of more than ", ROW_GROUP_ROWS, " rows.")
            part
        })
        ok <- TRUE
        for (column in definition$columns) {
            for (part in parts) {
                if (is.null(part[[column$name]])) {
                    if (isTRUE(column$required)) {
                        log$add("T0.1", where, "has no column `", column$name, "`.")
                        ok <- FALSE
                    }
                } else if (column_kind(part[[column$name]]) != column$type) {
                    log$add("T0.1", where, "The column `", column$name, "` is not of type ", column$type, ".")
                    ok <- FALSE
                } else if (!isTRUE(column$nullable) && column$type != "list<double>" && anyNA(part[[column$name]])) {
                    log$add("T0.1", where, "The column `", column$name, "` has missing values.")
                    ok <- FALSE
                }
            }
        }
        for (group in definition$one_of) {
            if (!any(unlist(group) %in% names(parts[[1]]))) log$add("T5.1", where, "needs one of the columns ", paste(unlist(group), collapse = ", "), ".")
        }
        for (column in definition$forbidden_columns) {
            if (column$name %in% names(parts[[1]])) log$add(column$rule, where, "must not have the column `", column$name, "`.")
        }
        # a column present under some gRNA integration strategies only (T3.3, T5.6)
        for (column in definition$columns) {
            if (is.null(column$present_under) || is.null(strategy)) next
            expected <- strategy %in% unlist(column$present_under)
            if (expected != (column$name %in% names(parts[[1]]))) {
                log$add(column$rule, where, if (expected) "needs" else "must not have", " the column `", column$name, "` under the strategy '", strategy, "'.")
            }
        }
        if (!ok) next
        for (part in parts) {
            factors <- vapply(part, is.factor, logical(1))
            part[factors] <- lapply(part[factors], as.character)
            key <- table_key(name, part)
            if (nrow(part) > 1L && !identical(byte_order(part, key), seq_len(nrow(part)))) log$add("T0.2", where, "The rows are not sorted by ", paste(key, collapse = ", "), ".")
        }
        whole <- read_table(tables[[name]])
        key <- table_key(name, whole)
        if (anyDuplicated(whole[key]) > 0L) log$add("T0.4", where, "Two rows have the same ", paste(key, collapse = ", "), ".")
        contents[[name]] <- whole
    }

    parameters <- records$parameters
    if (!is.null(parameters)) {
        if (!formula_in_subset(parameters$settings$formula)) log$add("T7.1", "stages/parameters.json", "The formula is outside the additive subset.")
        for (stage in c("calibration", "power", "discovery")) {
            fits <- contents[[paste0("fits/", stage)]]
            if (!is.null(fits) && any(lengths(fits$coefficients) != length(parameters$design$terms))) {
                log$add("T6.2", paste0("fits/", stage, ".parquet"), "A fit does not have one coefficient per term of the design.")
            }
        }
    }

    # against the data file
    if (!is.null(screen_file) && !is.null(analysis)) {
        screen <- screen_identity(screen_file, analysis$response_modality)
        if (is.null(screen$responses)) {
            log$add("A1.2", "analysis.json", "The data file has no modality named '", analysis$response_modality, "'.")
        } else {
            current <- list(
                id = screen$id, n_cells = length(screen$cells), n_responses = length(screen$responses), n_grnas = length(screen$grnas),
                cell_ids_sha256 = sha256_ids(screen$cells), response_ids_sha256 = sha256_ids(screen$responses), grna_ids_sha256 = sha256_ids(screen$grnas)
            )
            for (field in names(current)) {
                if (!isTRUE(all.equal(analysis$screen[[field]], current[[field]]))) log$add("A1.2", "analysis.json", "`screen.", field, "` is not that of the data file.")
            }
            known <- list(cell_id = screen$cells, response_id = screen$responses)
            for (name in names(contents)) {
                for (column in intersect(names(known), names(contents[[name]]))) {
                    if (!all(contents[[name]][[column]] %in% known[[column]])) log$add("T0.5", paste0(name, ".parquet"), "`", column, "` names something the data file does not have.")
                }
            }
            for (name in intersect(c("assignments", "negative_control_groups"), names(contents))) {
                if (!all(contents[[name]]$grna_id %in% screen$grnas)) log$add("T0.5", paste0(name, ".parquet"), "`grna_id` names something the data file does not have.")
            }
            for (name in intersect(c("pairs/discovery", "pairs/positive_control"), names(contents))) {
                if (!all(contents[[name]]$grna_target %in% screen$targets)) log$add("T0.5", paste0(name, ".parquet"), "`grna_target` names something the data file does not have.")
            }
            if (!is.null(contents$cells) && !setequal(contents$cells$cell_id, screen$cells)) log$add("T2.1", "cells.parquet", "does not have one row per cell of the data file.")
        }
    }
    log$all()
}
