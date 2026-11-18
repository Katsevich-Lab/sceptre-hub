# The sceptre file format (draft, format 0.1): R reader and writer.
#
# SPEC.md states the rules and RATIONALE.md the reasons. Section numbers below
# (S4, T3, ...) refer to both.
#
# Two layers:
#   1. Files <-> standard R containers, with no use of the sceptre package:
#        fmt_write_screen(), fmt_read_screen()   data file <-> SingleCellExperiment
#        fmt_init_analysis(), fmt_write_stage(), fmt_read_analysis(), fmt_clear_stage()
#   2. Checking files against SPEC.md, without building any object:
#        validate_screen(), validate_analysis()
#
# Nothing in this file uses the sceptre package. Building a sceptre_object from
# these files, and writing one, is in package-suite/sceptre_objects.R.
#
# The tables of an analysis, the stages and the shape of the JSON records are
# not written here: they are read from spec/ (tables.json, screen.json and two
# JSON Schemas), which the Python library reads too.
#
# Names. The data file names a modality for what was measured (`rna`), and the
# gRNA modality is always `grna` (S2). An analysis says which modality is its
# response (A1).

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
# The optional entries of a data file's metadata (S10.5): a name for people, and
# what is known about the screen at import and cannot be worked out later.
OPTIONAL_METADATA <- names(Filter(function(entry) !isTRUE(entry$required), SCREEN_SPEC$metadata$entries))
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

# The one type of both index arrays of a count matrix (S4.2, S4.3): 32-bit while
# the stored values, the cells and the features each number at most
# 2,147,483,647, and 64-bit beyond. It is the type scipy chooses for a matrix of
# that size. `index_limit` replaces the limit, for tests only: no test file is
# large enough to reach the real one.
index_type <- function(n_stored, n_cells, n_features, index_limit = NULL) {
    rule <- SCREEN_SPEC$counts$index
    limit <- if (is.null(index_limit)) rule$narrow_up_to else index_limit
    if (max(as.numeric(n_stored), as.numeric(n_cells), as.numeric(n_features)) <= limit) rule$narrow else rule$wide
}

# Counts as 32-bit integers and both index arrays in the one type the size
# requires, uncompressed (S4, S11). MuData would write doubles, so
# fmt_write_screen() has it write an empty matrix and puts
# the three arrays in afterwards. Rewriting full-size arrays in place would
# leave their first copy as dead space in the file: HDF5 does not give it back.
# `counts` is row-compressed features x cells, which is column-compressed
# cells x features in the file.
write_counts <- function(file, modality, counts, index_limit = NULL) {
    if (any(counts@x != round(counts@x))) fmt_stop("Counts must be integers.")
    group <- paste0("mod/", modality, "/X/")
    put <- function(name, values, h5type) {
        rhdf5::h5delete(file, paste0(group, name))
        rhdf5::h5createDataset(file, paste0(group, name), dims = length(values), H5type = h5type, chunk = NULL, level = 0)
        rhdf5::h5write(values, file, paste0(group, name))
    }
    put("data", as.integer(counts@x), "H5T_STD_I32LE")
    index <- c(int32 = "H5T_STD_I32LE", int64 = "H5T_STD_I64LE")[[index_type(length(counts@x), ncol(counts), nrow(counts), index_limit)]]
    # R has no 64-bit integer: doubles hold these whole numbers exactly, and HDF5 stores them in the type asked for
    as_stored <- if (index == "H5T_STD_I64LE") as.numeric else as.integer
    put("indices", as_stored(counts@j), index)
    put("indptr", as_stored(counts@p), index)
}

# Several gRNAs may be delivered on one vector, named by the optional column
# `vector_id` of the gRNA table (S8.6). A vector has one target, and a vector is
# not named like a gRNA on another vector. Returns the rules broken.
vector_faults <- function(grna_ids, vectors, targets) {
    faults <- list()
    if (any(tapply(targets, vectors, function(t) length(unique(t))) > 1L)) {
        faults <- c(faults, list(list(rule = "S8.7", message = "gRNAs on one vector have different targets; screens that combine targets on a vector are not supported.")))
    }
    own <- vectors[match(vectors, grna_ids)]   # the vector of the gRNA each vector is named like, if any
    if (any(!is.na(own) & own != vectors)) {
        faults <- c(faults, list(list(rule = "S8.8", message = "A vector is named like a gRNA on another vector.")))
    }
    faults
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
#                           `id` is kept if present; `name` (a label for people)
#                           and `source` (a URL, DOI or accession) are optional
fmt_write_screen <- function(sce, file, index_limit = NULL) {
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
    if (!is.null(grna_table$vector_id)) {
        if (anyNA(grna_table$vector_id)) fmt_stop("Every gRNA needs a `vector_id` when the column is present.")
        for (fault in vector_faults(rownames(grna), as.character(grna_table$vector_id), as.character(grna_table$grna_target))) fmt_stop(fault$message)
    }
    covariates <- split_covariates(SummarizedExperiment::colData(sce), c(main, GRNA))
    modality <- function(experiment, obs, var) {
        # every modality's cell table has its total count and its number of
        # features with a nonzero count; a writer computes them if they are absent (S6)
        counts <- SummarizedExperiment::assay(experiment, "counts")
        obs <- as.data.frame(obs, optional = TRUE, stringsAsFactors = FALSE)
        if (is.null(obs$n_umis)) obs$n_umis <- as.double(Matrix::colSums(counts))
        if (is.null(obs$n_nonzero)) obs$n_nonzero <- as.integer(Matrix::colSums(counts != 0))
        # row-compressed (features x cells) in R is column-compressed
        # (cells x features) in the file: one feature is one contiguous read (S3)
        stored <- methods::as(Matrix::drop0(methods::as(counts, "generalMatrix")), "RsparseMatrix")
        empty <- methods::as(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = dim(counts)), "RsparseMatrix")
        list(
            counts = stored,
            # MuData writes the layout with an empty matrix; write_counts() adds the arrays
            experiment = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = empty),
                colData = S4Vectors::DataFrame(as_file_frame(obs, TRUE), check.names = FALSE, row.names = cells),
                rowData = S4Vectors::DataFrame(as_file_frame(var, FALSE), check.names = FALSE, row.names = rownames(experiment))
            )
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
            modality(other, as.data.frame(SummarizedExperiment::colData(other), optional = TRUE), SummarizedExperiment::rowData(other))
        })
    )
    MuData::writeH5MU(
        MultiAssayExperiment::MultiAssayExperiment(
            lapply(experiments, function(e) e$experiment),
            colData = S4Vectors::DataFrame(as_file_frame(covariates$shared, TRUE), check.names = FALSE, row.names = cells),
            # a native record with scalar entries (S10)
            metadata = list(sceptre = c(
                list(
                    format = scalar(FORMAT_VERSION), moi = scalar(stamp$moi),
                    id = scalar(if (is.null(stamp$id)) new_id() else stamp$id)
                ),
                # optional: a label for people, and where the data came from (S10.5)
                lapply(stamp[intersect(OPTIONAL_METADATA, names(stamp)[!vapply(stamp, is.null, logical(1))])], scalar)
            ))
        ),
        file
    )
    for (name in names(experiments)) write_counts(file, name, experiments[[name]]$counts, index_limit)
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
    grna <- experiment(GRNA)
    # a reader refuses a file whose vectors break S8.7 or S8.8 (S8.9)
    grna_table <- SummarizedExperiment::rowData(grna)
    if (!is.null(grna_table$vector_id)) {
        for (fault in vector_faults(rownames(grna), as.character(grna_table$vector_id), as.character(grna_table$grna_target))) fmt_stop(file, ": ", fault$message)
    }
    SingleCellExperiment::altExp(sce, "gRNA") <- grna
    # modalities no analysis here uses are carried along, so that writing the
    # object back keeps them (S2.6)
    for (name in setdiff(names(mae), c(main, GRNA))) {
        other <- experiment(name)
        SummarizedExperiment::colData(other) <- S4Vectors::DataFrame(
            cell_table(SummarizedExperiment::colData(mae[[name]])), check.names = FALSE, row.names = cells
        )
        SingleCellExperiment::altExp(sce, name) <- other
    }
    S4Vectors::metadata(sce) <- list(sceptre = c(
        list(format = as.character(stamp$format), moi = stamp$moi, id = stamp$id, response_modality = main),
        stamp[intersect(OPTIONAL_METADATA, names(stamp))]
    ))
    sce
}

# What an analysis records about its data file (A1).
screen_link <- function(sce, file = NULL, directory = NULL) {
    grna <- SingleCellExperiment::altExp(sce, "gRNA")
    link <- list(
        id = S4Vectors::metadata(sce)$sceptre$id,
        n_cells = ncol(sce), n_responses = nrow(sce), n_grnas = nrow(grna),
        cell_ids_sha256 = sha256_ids(colnames(sce)),
        response_ids_sha256 = sha256_ids(rownames(sce)),
        grna_ids_sha256 = sha256_ids(rownames(grna))
    )
    # for people: which dataset this is; not part of the identity a reader checks
    link$name <- S4Vectors::metadata(sce)$sceptre$name
    # a hint only, and only when the data file is near the analysis (A1.4)
    if (!is.null(file)) link$path <- nearby_path(file, directory)
    link
}
# The data file's path relative to an analysis folder, if the file is in the
# folder's parent or grandparent (the `analyses/<name>/` convention) or below
# one of them; otherwise nothing. A far-away file would give a long chain of
# `..` that names directories and helps nobody.
nearby_path <- function(file, directory) {
    path <- relative_path(file, directory)
    if (sum(strsplit(path, "/", fixed = TRUE)[[1]] == "..") > 2L) NULL else path
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

# An analysis folder has one folder per stage: the stage's record and its
# tables (A2). A table's name is its path without `.parquet`: "qc/cells".
stage_file <- function(directory, stage) file.path(directory, stage, "record.json")
fmt_stages_run <- function(directory) STAGES[file.exists(stage_file(directory, STAGES))]
table_path <- function(directory, name) file.path(directory, paste0(name, ".parquet"))

# Write one stage into its folder: its tables first, then its record, which
# lists them (A2, A5). `record` holds settings, and optionally summaries.
fmt_write_stage <- function(directory, stage, record, tables = list(), provenance = NULL) {
    if (!file.exists(file.path(directory, "analysis.json"))) fmt_stop(directory, " is not an analysis folder.")
    earlier <- names(STAGE_RANK)[STAGE_RANK < STAGE_RANK[[stage]]]
    missing <- setdiff(earlier, fmt_stages_run(directory))
    if (length(missing) > 0L) fmt_stop("Stage '", stage, "' needs the stage(s) ", paste(missing, collapse = ", "), " first.")
    if (file.exists(stage_file(directory, stage))) fmt_stop("Stage '", stage, "' is already written; clear it first.")
    # a folder without a record is a stage that did not finish: nothing of it is kept (A4.1)
    unlink(file.path(directory, stage), recursive = TRUE)
    files <- lapply(names(tables), function(name) {
        if (is.null(TABLES[[name]]) || !identical(TABLES[[name]]$stage, stage)) fmt_stop("Stage '", stage, "' does not write the table ", name, ".")
        # listed by its name within the stage's folder
        list(path = paste0(basename(name), ".parquet"), rows = write_table(tables[[name]], table_path(directory, name), table_key(name, tables[[name]])))
    })
    if (is.null(provenance)) provenance <- list(package = "sceptre_format.R", version = "prototype", language = "R")
    if (is.null(provenance$written_at)) provenance$written_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    write_json(c(list(stage = stage), record, list(provenance = provenance, files = files)), stage_file(directory, stage))
    invisible(directory)
}

# Remove a stage and every later stage: their folders, records first (A4). The
# fit cache goes too when the QC stage or an earlier one is cleared, and stays
# when only a testing stage is (T6).
fmt_clear_stage <- function(directory, stage) {
    cleared <- names(STAGE_RANK)[STAGE_RANK > STAGE_RANK[[stage]] | names(STAGE_RANK) == stage]
    for (later in cleared) {
        # the record goes first: a stage without a record has not run
        unlink(stage_file(directory, later))
        unlink(file.path(directory, later), recursive = TRUE)
    }
    if (TABLE_SPEC$fits$cleared_with %in% cleared) unlink(file.path(directory, FITS), recursive = TRUE)
    invisible(directory)
}

# ---- the fit cache (T6): fits/response/, a folder of part files outside the
# stages. It is no stage's output and is listed in no record. Each writer adds
# a part of its own; what makes the fits valid is recorded once, beside them.
FITS <- TABLE_SPEC$fits$folder
fits_identity_file <- function(directory) file.path(directory, TABLE_SPEC$fits$identity)
# A formula is the same formula however it is spaced (T6.5).
formula_key <- function(text) gsub("[[:space:]]", "", as.character(text))
identity_value <- function(identity, name) {
    value <- as.character(identity[[name]])
    if (identical(name, "formula")) formula_key(value) else value
}
same_identity <- function(a, b) {
    all(vapply(union(names(a), names(b)), function(name) identical(identity_value(a, name), identity_value(b, name)), logical(1)))
}
fmt_write_fits <- function(directory, table, identity) {
    missing <- setdiff(c("parameters", "qc"), fmt_stages_run(directory))
    if (length(missing) > 0L) fmt_stop("Fits can be saved once the parameters and qc stages are written.")
    path <- fits_identity_file(directory)
    if (file.exists(path)) {
        # a later writer checks the record; it does not rewrite it
        if (!same_identity(read_json(path), identity)) fmt_stop("The fits already saved in ", FITS, " were made differently.")
    } else {
        write_json(identity, path)
    }
    part <- file.path(directory, FITS, paste0("part-", substr(new_id(), 1L, 12L), ".parquet"))
    write_table(table, part, table_key(FITS, table))
    invisible(part)
}
# The cached fits and how they were made, or NULL. Two parts may hold the same
# response; a reader takes either.
fmt_read_fits <- function(directory) {
    if (!file.exists(fits_identity_file(directory))) return(NULL)
    folder <- file.path(directory, FITS)
    table <- if (length(table_files(folder)) > 0L) read_table(folder, TABLES[[FITS]]$columns) else NULL
    if (!is.null(table)) table <- table[!duplicated(table$response_id), , drop = FALSE]
    list(identity = read_json(fits_identity_file(directory)), table = table)
}

# Read an analysis folder: its link, the record of each stage that ran, their
# tables (named as in TABLES), and the fit cache if there is one.
fmt_read_analysis <- function(directory) {
    analysis <- read_json(file.path(directory, "analysis.json"))
    check_format(analysis$format, directory)
    stages <- lapply(stats::setNames(nm = fmt_stages_run(directory)), function(stage) read_json(stage_file(directory, stage)))
    tables <- list()
    for (stage in names(stages)) {
        record <- stages[[stage]]
        # a record is the record of its folder's stage, and lists files of that folder only (A2.1, A5.3)
        if (!identical(record$stage, stage)) fmt_stop(stage, "/record.json says it is the record of the stage '", record$stage, "'.")
        for (file in record$files) {
            if (!grepl("^[A-Za-z0-9_]+\\.parquet$", file$path)) fmt_stop(stage, "/record.json lists '", file$path, "', which is not a table of its own folder.")
            name <- paste0(stage, "/", sub("\\.parquet$", "", file$path))
            path <- table_path(directory, name)
            if (!file.exists(path)) fmt_stop("The stage '", record$stage, "' lists ", file$path, ", which is missing.")
            if (table_rows(path) != file$rows) {
                fmt_stop(name, ".parquet has ", table_rows(path), " rows; the stage '", record$stage, "' wrote ", file$rows, ".")
            }
            tables[[name]] <- read_table(path, TABLES[[name]]$columns)
        }
    }
    fits <- fmt_read_fits(directory)
    tables[[FITS]] <- fits$table
    list(
        format = analysis$format, response_modality = analysis$response_modality,
        screen = analysis$screen, stages = stages, tables = tables, fits_identity = fits$identity
    )
}

# Refuse a data file other than the one the analysis was run on (A1).
check_link <- function(recorded, sce, what) {
    current <- screen_link(sce)
    for (field in c("id", "cell_ids_sha256", "response_ids_sha256", "grna_ids_sha256")) {
        if (!identical(recorded[[field]], current[[field]])) {
            named <- if (is.null(recorded$name)) "" else paste0(", the dataset '", recorded$name, "'")
            fmt_stop("The analysis ", what, " was run on another data file", named, " (its `", field, "` differs).")
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

# ---------------------------------------------------------------- comparing and describing

rows_sorted <- function(table) {
    table <- as.data.frame(table, stringsAsFactors = FALSE)
    if (nrow(table) == 0L) return(table)
    table <- table[do.call(order, c(unname(as.list(table)), list(method = "radix"))), , drop = FALSE]
    rownames(table) <- NULL
    table
}

# The files of a folder with their sizes.
file_listing <- function(directory) {
    files <- list.files(directory, recursive = TRUE)
    data.frame(file = files, kilobytes = round(file.size(file.path(directory, files)) / 1024, 1))
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

# `index_limit` is for tests only: see index_type().
validate_screen <- function(file, index_limit = NULL) {
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
    for (name in names(entries)[!vapply(entries, function(e) isTRUE(e$required), logical(1))]) {
        if (nrow(entry(name)) == 1L && !scalar_of(name, classes[[name]])) log$add("S10.5", paste0(where, "/", name), "must be a scalar text entry.")
    }
    required <- vapply(entries, function(e) isTRUE(e$required), logical(1))
    if (!all(mapply(scalar_of, names(entries)[required], classes[required]))) {
        log$add("S10.1", where, "needs the scalar text entries `format`, `moi` and `id`.")
    } else if (is.null(format_parts(read(paste0(where, "/format"))))) {
        log$add("S10.1", paste0(where, "/format"), "is not a version of the form MAJOR.MINOR.")
    } else {
        if (!read(paste0(where, "/moi")) %in% unlist(entries$moi$values)) log$add("S10.2", paste0(where, "/moi"), "must be 'high' or 'low'.")
        if (!known_major(read(paste0(where, "/format")))) log$add("S10.3", paste0(where, "/format"), "is a format this validator does not know.")
        if (!nzchar(read(paste0(where, "/id")))) log$add("S10.1", paste0(where, "/id"), "is empty.")
    }
    # optional entries with a fixed vocabulary (S10.6)
    for (name in OPTIONAL_METADATA) {
        allowed <- unlist(entries[[name]]$values)
        if (!is.null(allowed) && nrow(entry(name)) == 1L && scalar_of(name, "STRING") && !read(paste0(where, "/", name)) %in% allowed) {
            log$add("S10.6", paste0(where, "/", name), "must be one of ", paste(allowed, collapse = ", "), ".")
        }
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
        # both index arrays in the one type this file's size requires (S4.2, S4.3)
        index <- index_type(as.numeric(stored$dim[stored$name == "data"]), length(cells), n_features, index_limit)
        for (array in names(expect)) {
            row <- stored[stored$name == array, , drop = FALSE]
            kind <- if (expect[[array]] == "index") index else expect[[array]]
            why <- if (expect[[array]] == "index") " (the index type a file of this size requires)" else ""
            if (!is_stored_as(row$dtype, kind)) {
                # another whole-number type: not conforming, and a reader still accepts it (S4.4)
                add <- if (row$dclass %in% c("INTEGER", "FLOAT")) log$add_readable else log$add
                add(rules[[array]], paste0(x, "/", array), "is not stored as ", sub("int", "", kind), "-bit signed integers", why, ".")
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
            } else if (isTRUE(column$required) || isTRUE(column$no_missing)) {
                values <- text(column$name)
                if (anyNA(values) || any(!nzchar(values))) log$add(column$rule, path, "must be text with no missing values.")
            }
        }
        # target positions: 1-based, both ends included; all three or none (S8.10, S8.11)
        position <- function(name, as_text) {
            row <- columns[columns$name == name, , drop = FALSE]
            if (nrow(row) == 0L) return(rep(NA, length(features[[GRNA]])))
            kind <- storage_kind(as.vector(attribute(paste0(table, "/", name), "encoding-type")), row$dtype)
            if (as_text) {
                if (kind %in% c("string", "categorical")) text(name) else NULL
            } else {
                if (kind %in% c("double", "int32")) suppressMessages(as.numeric(read(paste0(table, "/", name)))) else NULL
            }
        }
        chromosome <- position("target_chr", TRUE)
        start <- position("target_start", FALSE)
        end <- position("target_end", FALSE)
        if (!is.null(chromosome) && !is.null(start) && !is.null(end)) {
            given <- cbind(!is.na(chromosome), !is.na(start), !is.na(end))
            if (any(rowSums(given) %in% 1:2)) log$add("S8.11", table, "A gRNA has some of `target_chr`, `target_start`, `target_end` and not all three.")
            both <- !is.na(start) & !is.na(end)
            if (any(start[both] != round(start[both]) | end[both] != round(end[both]) | start[both] < 1 | end[both] < start[both])) {
                log$add("S8.10", table, "Target positions must be whole numbers, 1-based with both ends included: `target_start` at least 1 and `target_end` not before it.")
            }
        }
        # gRNAs delivered on one vector (S8.6 to S8.8)
        if ("vector_id" %in% columns$name && "grna_target" %in% columns$name) {
            vectors <- tryCatch(text("vector_id"), error = function(e) NULL)
            targets <- tryCatch(text("grna_target"), error = function(e) NULL)
            if (!is.null(vectors) && !is.null(targets) && !anyNA(vectors) && !anyNA(targets)) {
                for (fault in vector_faults(features[[GRNA]], vectors, targets)) log$add(fault$rule, paste0(table, "/vector_id"), fault$message)
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
    # the form of the hint has a rule of its own
    for (fault in faults) log$add(if (startsWith(fault, "/screen/path ")) "A1.4" else "A1.1", "analysis.json", fault)
    analysis <- if (length(faults) == 0L) read_json(path) else NULL
    if (!is.null(analysis) && !known_major(analysis$format)) log$add("A1.1", "analysis.json", "is in a format this validator does not know.")

    # stage records: <stage>/record.json
    present <- fmt_stages_run(directory)
    records <- list()
    unchecked <- character()   # stages whose record cannot be trusted
    for (stage in present) {
        where <- paste0(stage, "/record.json")
        faults <- schema_problems(stage_file(directory, stage), "stage.schema.json")
        for (fault in faults) log$add("A2.1", where, fault)
        if (length(faults) > 0L) {
            unchecked <- c(unchecked, stage)
            next
        }
        record <- read_json(stage_file(directory, stage))
        if (!identical(record$stage, stage)) {
            log$add("A2.1", where, "says it is the record of the stage '", record$stage, "'.")
            unchecked <- c(unchecked, stage)
            next
        }
        records[[stage]] <- record
    }
    for (stage in present) {
        missing <- setdiff(names(STAGE_RANK)[STAGE_RANK < STAGE_RANK[[stage]]], present)
        if (length(missing) > 0L) log$add("A2.2", paste0(stage, "/record.json"), "The stage needs the stage(s) ", paste(missing, collapse = ", "), " first.")
    }

    # the files each record lists, all in the stage's own folder
    listed <- character()
    tables <- list()
    for (stage in names(records)) {
        record <- records[[stage]]
        where <- paste0(stage, "/record.json")
        names_listed <- paste0(stage, "/", sub("\\.parquet$", "", vapply(record$files, function(f) f$path, character(1))))
        listed <- c(listed, names_listed)
        for (file in record$files) {
            name <- paste0(stage, "/", sub("\\.parquet$", "", file$path))
            shown <- paste0(name, ".parquet")
            full <- table_path(directory, name)
            if (!file.exists(full)) {
                log$add("A5.1", shown, "is listed by its stage's record and is missing.")
                next
            }
            rows <- tryCatch(table_rows(full), error = function(e) NA)
            if (is.na(rows)) {
                log$add("A5.1", shown, "cannot be read as Parquet.")
                next
            }
            if (rows != file$rows) log$add("A5.2", shown, "has ", rows, " rows; the stage '", stage, "' wrote ", file$rows, ".")
            if (!is.null(TABLES[[name]])) tables[[name]] <- full
        }
        required <- names(TABLES)[vapply(TABLES, function(t) identical(t$stage, stage) && isTRUE(t$required), logical(1))]
        for (name in setdiff(required, names_listed)) log$add("A5.4", where, "does not list the table ", name, ".")
    }
    # a stage's table that is in its folder and that no record lists (A4.1)
    for (name in names(TABLES)) {
        stage <- TABLES[[name]]$stage
        if (is.null(stage) || stage %in% unchecked || name %in% listed) next
        if (file.exists(table_path(directory, name))) log$add("A4.1", paste0(name, ".parquet"), "is in the folder and its stage's record does not list it.")
    }

    # the fit cache: parts, and one record of how the fits were made (T6)
    fits_identity <- NULL
    fits_folder <- file.path(directory, FITS)
    if (dir.exists(fits_folder)) {
        if (!file.exists(fits_identity_file(directory))) {
            log$add("T6.1", FITS, "has no record of how the fits were made.")
        } else {
            faults <- schema_problems(fits_identity_file(directory), "fits.schema.json")
            for (fault in faults) log$add("T6.1", TABLE_SPEC$fits$identity, fault)
            if (length(faults) == 0L) fits_identity <- read_json(fits_identity_file(directory))
        }
        if (length(table_files(fits_folder)) > 0L) tables[[FITS]] <- fits_folder
        if (!all(c("parameters", "qc") %in% present)) log$add("T6.3", FITS, "The fit cache needs the parameters and qc stages.")
    }

    # each table against its definition
    strategy <- records$parameters$settings$grna_integration_strategy
    contents <- list()
    for (name in names(tables)) {
        definition <- TABLES[[name]]$definition
        where <- if (identical(name, FITS)) name else paste0(name, ".parquet")
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
        # within one file the key is sorted and unique; across the parts of the
        # fit cache it need not be unique
        across_parts <- !identical(definition$unique_across_parts, FALSE)
        for (part in parts) {
            factors <- vapply(part, is.factor, logical(1))
            part[factors] <- lapply(part[factors], as.character)
            key <- table_key(name, part)
            if (nrow(part) > 1L && !identical(byte_order(part, key), seq_len(nrow(part)))) log$add("T0.2", where, "The rows are not sorted by ", paste(key, collapse = ", "), ".")
            if (!across_parts && anyDuplicated(part[key]) > 0L) log$add("T0.4", where, "Two rows of one file have the same ", paste(key, collapse = ", "), ".")
        }
        # the parts already read, put together: each table is read once
        plain <- lapply(parts, function(part) {
            factors <- vapply(part, is.factor, logical(1))
            part[factors] <- lapply(part[factors], as.character)
            part
        })
        shared <- Reduce(intersect, lapply(plain, names))
        whole <- do.call(rbind, lapply(plain, function(part) part[shared]))
        rownames(whole) <- NULL
        key <- table_key(name, whole)
        if (across_parts && anyDuplicated(whole[key]) > 0L) log$add("T0.4", where, "Two rows have the same ", paste(key, collapse = ", "), ".")
        contents[[name]] <- whole
    }

    parameters <- records$parameters
    if (!is.null(parameters)) {
        if (!formula_in_subset(parameters$settings$formula)) log$add("T7.1", "parameters/record.json", "The formula is outside the additive subset.")
        fits <- contents[[FITS]]
        if (!is.null(fits) && any(lengths(fits$coefficients) != length(parameters$design$terms))) {
            log$add("T6.2", FITS, "A fit does not have one coefficient per term of the design.")
        }
        # the cache was made for this analysis's parameters (T6.3)
        if (!is.null(fits_identity)) {
            expected <- list(
                formula = parameters$settings$formula, response_fit_method = parameters$settings$response_fit_method,
                control_group = parameters$settings$control_group, screen_id = analysis$screen$id
            )
            for (field in names(expected)) {
                if (!is.null(expected[[field]]) && !identical(identity_value(fits_identity, field), identity_value(expected, field))) {
                    log$add("T6.3", TABLE_SPEC$fits$identity, "`", field, "` is not that of the analysis.")
                }
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
            for (name in intersect(c("assignment/assignments", "calibration/negative_control_groups"), names(contents))) {
                if (!all(contents[[name]]$grna_id %in% screen$grnas)) log$add("T0.5", paste0(name, ".parquet"), "`grna_id` names something the data file does not have.")
            }
            for (name in intersect(c("parameters/discovery_pairs", "parameters/positive_control_pairs"), names(contents))) {
                if (!all(contents[[name]]$grna_target %in% screen$targets)) log$add("T0.5", paste0(name, ".parquet"), "`grna_target` names something the data file does not have.")
            }
            kept <- contents[["qc/cells"]]
            if (!is.null(kept) && !setequal(kept$cell_id, screen$cells)) {
                log$add("T2.1", "qc/cells.parquet", "does not have one row per cell of the data file.")
            } else if (!is.null(kept) && !is.null(fits_identity)) {
                in_use <- screen$cells[screen$cells %in% kept$cell_id[kept$in_use]]
                if (!identical(fits_identity$cells_in_use_sha256, sha256_ids(in_use))) log$add("T6.3", TABLE_SPEC$fits$identity, "`cells_in_use_sha256` is not that of the cells in use.")
            }
        }
    }
    log$all()
}
