# The sceptre analysis format (draft, format 0.1): R validator, and temporary
# writers and readers.
#
# SPEC.md states the rules and RATIONALE.md the reasons. Section numbers below
# (A1, T3, ...) refer to both. The tables of an analysis, the stages and the
# shape of the JSON records are read from spec/ (tables.json and three JSON
# Schemas), which the Python code reads too.
#
# This builds on the Perturb-seq data file format, whose code format/setup.R
# sources first: an analysis names its data file, and checking an analysis
# fully needs that file.
#
#   fmt_init_analysis(), fmt_write_stage(), fmt_read_analysis(), fmt_clear_stage()
#   validate_analysis()   checks an analysis folder against SPEC.md
#
# Nothing here uses the sceptre package.

# The version of the analysis format, "MAJOR.MINOR" (A1.1). It is separate from
# the data file format's version.
ANALYSIS_FORMAT_VERSION <- "0.1"
ANALYSIS_SPEC_DIR <- file.path(Sys.getenv("SCEPTRE_ANALYSIS_FORMAT_DIR", getwd()), "spec")

TABLE_SPEC <- jsonlite::fromJSON(file.path(ANALYSIS_SPEC_DIR, "tables.json"), simplifyVector = FALSE)

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
# `metadata` is text by key, stored in the file's own key-value metadata.
write_table <- function(table, path, key, metadata = NULL) {
    table <- as.data.frame(table, stringsAsFactors = FALSE)
    factors <- vapply(table, is.factor, logical(1))
    table[factors] <- lapply(table[factors], as.character)
    table <- table[byte_order(table, key), , drop = FALSE]
    rownames(table) <- NULL
    atomically(path, function(temporary) {
        nanoparquet::write_parquet(
            table, temporary, metadata = metadata,
            options = nanoparquet::parquet_options(num_rows_per_row_group = ROW_GROUP_ROWS)
        )
    })
    nrow(table)
}

table_files <- function(path) {
    if (length(path) == 1L && dir.exists(path)) sort(list.files(path, pattern = "\\.parquet$", full.names = TRUE)) else path
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

# What an analysis records about its data file (A1).
screen_link <- function(sce, file = NULL, directory = NULL) {
    grna <- SingleCellExperiment::altExp(sce, "gRNA")
    link <- list(
        id = S4Vectors::metadata(sce)$dataset_id,
        # the data file format the analysis was written against; recorded, not checked
        format = S4Vectors::metadata(sce)$perturb_seq_format,
        n_cells = ncol(sce), n_responses = nrow(sce), n_grnas = nrow(grna),
        cell_ids_sha256 = sha256_ids(colnames(sce)),
        response_ids_sha256 = sha256_ids(rownames(sce)),
        grna_ids_sha256 = sha256_ids(rownames(grna))
    )
    # for people: which dataset this is; not part of the identity a reader checks
    link$title <- S4Vectors::metadata(sce)$title
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
        list(format = ANALYSIS_FORMAT_VERSION, response_modality = response_modality, screen = link),
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
    if (is.null(provenance)) provenance <- list(package = "sceptre_analysis.R", version = "prototype", language = "R")
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

# ---- the fit cache (T6): fits/response/, a folder of part files and nothing
# else, outside the stages. It is no stage's output and is listed in no record.
# Each writer adds a part of its own, and each part says how its fits were made,
# in the Parquet file's own key-value metadata. No file is shared between
# writers, so there is nothing to keep in step: a part made for other cells,
# another formula or another method is simply not used.
FITS <- TABLE_SPEC$fits$folder

FIT_KEY <- TABLE_SPEC$fits$identity_key

# A formula is the same formula however it is spaced (T6.5).
formula_key <- function(text) gsub("[[:space:]]", "", as.character(text))

identity_value <- function(identity, name) {
    value <- as.character(identity[[name]])
    if (identical(name, "formula")) formula_key(value) else value
}

# Identities are compared as parsed values, never as text.
same_identity <- function(a, b) {
    all(vapply(union(names(a), names(b)), function(name) identical(identity_value(a, name), identity_value(b, name)), logical(1)))
}

# The identity a part carries, as the text stored, or NULL.
fit_identity_text <- function(file) {
    stored <- nanoparquet::read_parquet_metadata(file)$file_meta_data$key_value_metadata[[1]]
    value <- stored$value[stored$key == FIT_KEY]
    if (length(value) == 1L && !is.na(value)) value else NULL
}

# Every part of the cache with the identity it carries. The identity is NULL
# where a part has none or it is not a JSON object; a reader ignores such a
# part (T6.5).
fit_parts <- function(directory) {
    folder <- file.path(directory, FITS)
    files <- if (dir.exists(folder)) table_files(folder) else character()
    lapply(files, function(file) {
        text <- fit_identity_text(file)
        identity <- if (is.null(text)) NULL else tryCatch(jsonlite::fromJSON(text, simplifyVector = TRUE, simplifyDataFrame = FALSE), error = function(e) NULL)
        list(file = file, identity = if (is.list(identity) && !is.null(names(identity))) identity else NULL)
    })
}

# What the folder alone says of the identity its fits must have: everything but
# the cells in use, whose hash needs the data file's order of cells.
analysis_fit_fields <- function(directory) {
    if (!all(c("parameters", "qc") %in% fmt_stages_run(directory))) return(NULL)
    settings <- read_json(stage_file(directory, "parameters"))$settings
    list(
        screen_id = read_json(file.path(directory, "analysis.json"))$screen$id, formula = settings$formula,
        response_fit_method = settings$response_fit_method, control_group = settings$control_group
    )
}

identity_matches <- function(identity, fields) {
    !is.null(identity) && all(vapply(names(fields), function(name) identical(identity_value(identity, name), identity_value(fields, name)), logical(1)))
}

# Add a part holding `table`, with the identity of its fits inside it.
fmt_write_fits <- function(directory, table, identity) {
    missing <- setdiff(c("parameters", "qc"), fmt_stages_run(directory))
    if (length(missing) > 0L) fmt_stop("Fits can be saved once the parameters and qc stages are written.")
    part <- file.path(directory, FITS, paste0("part-", substr(new_id(), 1L, 12L), ".parquet"))
    text <- as.character(jsonlite::toJSON(identity[sort(names(identity))], auto_unbox = TRUE))
    write_table(table, part, table_key(FITS, table), stats::setNames(text, FIT_KEY))
    invisible(part)
}

# The cached fits made as `identity` says, and that identity, or NULL if no
# part has it. Parts with another identity, or with none, are ignored. Two parts
# may hold the same response; a reader takes either.
# Without `identity`, the parts made for this analysis's own parameters are
# returned. If those were made for more than one set of cells (parts left by an
# earlier QC), the folder alone cannot say which set is current, and NULL is
# returned: give the identity.
fmt_read_fits <- function(directory, identity = NULL) {
    parts <- fit_parts(directory)
    if (is.null(identity)) {
        fields <- analysis_fit_fields(directory)
        if (is.null(fields)) return(NULL)
        found <- Filter(function(part) identity_matches(part$identity, fields), parts)
        if (length(found) == 0L || !all(vapply(found, function(part) same_identity(part$identity, found[[1]]$identity), logical(1)))) return(NULL)
        identity <- found[[1]]$identity
    }
    files <- unlist(lapply(parts, function(part) if (!is.null(part$identity) && same_identity(part$identity, identity)) part$file))
    if (length(files) == 0L) return(NULL)
    table <- read_table(files, TABLES[[FITS]]$columns)
    list(identity = identity, table = table[!duplicated(table$response_id), , drop = FALSE])
}

# Read an analysis folder: its link, the record of each stage that ran, their
# tables (named as in TABLES), and the fit cache if there is one.
fmt_read_analysis <- function(directory) {
    analysis <- read_json(file.path(directory, "analysis.json"))
    check_format(analysis$format, directory, ANALYSIS_FORMAT_VERSION)
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
            named <- if (is.null(recorded$title)) "" else paste0(", the dataset '", recorded$title, "'")
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

# What validate_analysis() needs from a data file, read without building an object.
screen_identity <- function(file, response_modality) {
    read <- function(path) as.character(as.vector(rhdf5::h5read(file, path)))
    index <- function(table) read(paste0(table, "/", as.vector(rhdf5::h5readAttributes(file, table)[["_index"]])))
    modalities <- subset(rhdf5::h5ls(file, recursive = 2L), group == "/mod")$name
    out <- list(modalities = modalities, id = read("uns/dataset_id"), cells = index("mod/grna/obs"), grnas = index("mod/grna/var"))
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
        if (is.null(cache[[name]])) cache[[name]] <<- jsonvalidate::json_validator(file.path(ANALYSIS_SPEC_DIR, name), engine = "ajv")
        cache[[name]]
    }
})

# Problems a JSON file has against a schema, as text; character() when it conforms.
schema_problems <- function(path, schema) schema_problems_text(paste(readLines(path, warn = FALSE), collapse = "\n"), schema)

schema_problems_text <- function(text, schema) {
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
    if (!is.null(analysis) && !known_major(analysis$format, ANALYSIS_FORMAT_VERSION)) log$add("A1.1", "analysis.json", "is in a format this validator does not know.")

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

    # the fit cache: part files only, each saying how its fits were made (T6)
    fit_identities <- list()
    fits_folder <- file.path(directory, FITS)
    if (dir.exists(fits_folder)) {
        entries <- sort(list.files(fits_folder, all.files = TRUE, no.. = TRUE), method = "radix")
        for (entry in entries) {
            if (dir.exists(file.path(fits_folder, entry)) || !grepl("\\.parquet$", entry)) log$add("T6.1", paste0(FITS, "/", entry), "is not a part file; the fit cache holds part files only.")
        }
        for (file in table_files(fits_folder)) {
            text <- tryCatch(list(fit_identity_text(file)), error = function(e) NULL)
            if (is.null(text)) next   # not a Parquet file: reported with the table below
            where <- paste0(FITS, "/", basename(file))
            if (is.null(text[[1]])) {
                log$add("T6.1", where, "does not say how its fits were made: its metadata has no `", FIT_KEY, "`.")
                next
            }
            faults <- schema_problems_text(text[[1]], TABLE_SPEC$fits$identity_schema)
            for (fault in faults) log$add("T6.1", where, "`", FIT_KEY, "` ", fault)
            if (length(faults) == 0L) fit_identities[[file]] <- jsonlite::fromJSON(text[[1]], simplifyVector = TRUE, simplifyDataFrame = FALSE)
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
        # a stored logarithm agrees with the column it is the logarithm of (T5.2)
        agree <- definition$log_2_of
        if (!is.null(agree) && ok && agree$of %in% names(parts[[1]])) {
            for (part in parts) {
                held <- part[[agree$column]]
                derived <- log_2_of(part[[agree$of]])
                both <- !is.na(held) & !is.na(derived)
                if (any(held[both] != derived[both] & abs(held[both] - derived[both]) > agree$relative_tolerance * abs(derived[both]))) {
                    log$add(agree$rule, where, "`", agree$column, "` is not the base-2 logarithm of `", agree$of, "`.")
                    break
                }
            }
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
        # A part made for this analysis's parameters has one coefficient per term of its
        # design (T6.2). A part made otherwise is allowed and unused, and is not
        # reported: nothing depends on it, and clearing it is housekeeping (T6.3, T6.6).
        expected <- list(
            formula = parameters$settings$formula, response_fit_method = parameters$settings$response_fit_method,
            control_group = parameters$settings$control_group
        )
        if (!is.null(analysis)) expected$screen_id <- analysis$screen$id
        for (file in names(fit_identities)) {
            if (!identity_matches(fit_identities[[file]], expected)) next
            part <- as.data.frame(nanoparquet::read_parquet(file), stringsAsFactors = FALSE)
            if (!is.null(part$coefficients) && any(lengths(part$coefficients) != length(parameters$design$terms))) {
                log$add("T6.2", FITS, "A fit does not have one coefficient per term of the design.")
                break
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
            }
        }
    }
    log$all()
}
