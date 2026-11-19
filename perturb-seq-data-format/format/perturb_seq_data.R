# The Perturb-seq data file format (draft, format 0.1): R validator, and a
# temporary writer.
#
# SPEC.md states the rules and RATIONALE.md the reasons. Section numbers below
# (S4, S10, ...) refer to both. The structure of a data file is read from
# spec/screen.json, which the Python code reads too.
#
#   validate_screen()    checks a data file against SPEC.md, building no object
#   check_screen()       says whether software must refuse a file it has read
#   fmt_write_screen()   SingleCellExperiment -> data file; MuData::readH5MU
#                        reads it back
#
# Nothing here knows about analyses or about any analysis tool.
#
# Names. A modality is named for what was measured (`rna`), and the gRNA
# modality is always `grna` (S2).

# The format version, "MAJOR.MINOR" (S10.3). A minor increase adds things an
# older reader can ignore; a major increase is one an older reader must refuse.
FORMAT_VERSION <- "0.1"

format_parts <- function(version) {
    version <- as.character(as.vector(version))
    if (length(version) != 1L || is.na(version) || !grepl("^[0-9]+\\.[0-9]+$", version)) return(NULL)
    as.integer(strsplit(version, ".", fixed = TRUE)[[1]])
}

known_major <- function(version, known = FORMAT_VERSION) identical(format_parts(version)[1], format_parts(known)[1])

SPEC_DIR <- file.path(Sys.getenv("PERTURB_SEQ_DATA_FORMAT_DIR", getwd()), "spec")

SCREEN_SPEC <- jsonlite::fromJSON(file.path(SPEC_DIR, "screen.json"), simplifyVector = FALSE)

# The optional entries of a data file's metadata (S10.5): a name for people, and
# what is known about the screen at import and cannot be worked out later.
OPTIONAL_METADATA <- names(Filter(function(entry) !isTRUE(entry$required), SCREEN_SPEC$metadata$entries))

GRNA <- "grna"

# the name fmt_write_screen() gives the main experiment's modality when the
# object does not name it
DEFAULT_MAIN_MODALITY <- "rna"

# the columns of the gRNA table that hold a target's position (S8.2)
GRNA_COORDINATES <- c(chr = "target_chr", start = "target_start", end = "target_end")

fmt_stop <- function(...) stop(..., call. = FALSE)

# ---------------------------------------------------------------- small helpers

new_id <- function() unname(tools::md5sum(bytes = charToRaw(paste(Sys.time(), Sys.getpid(), tempfile(), stats::runif(1)))))

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

# The writer's input has one table of cell covariates, a modality's columns
# prefixed with `<modality>_`. It is split in the file (S5): facts about the
# cell are in the shared cell table, measurements of a modality in that
# modality's, without the prefix.
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

# The optional column `grna_type` of the gRNA table says what kind of gRNA each
# is (S8.15), and agrees with `grna_target`: a control of any kind has the target
# `non-targeting`, and a targeting gRNA has another target (S8.16). A control may
# have a position: where it cuts. Returns the rules broken.
GRNA_TYPES <- unlist(Filter(function(column) column$name == "grna_type", SCREEN_SPEC$grna_table)[[1]]$values)
type_faults <- function(types, targets) {
    if (anyNA(types) || !all(types %in% GRNA_TYPES)) {
        return(list(list(rule = "S8.15", message = paste0("`grna_type` must be one of ", paste(GRNA_TYPES, collapse = ", "), "."))))
    }
    control <- types != "targeting"
    if (any(control != (targets == "non-targeting"))) {
        return(list(list(rule = "S8.16", message = "`grna_type` and `grna_target` disagree: a control of any kind has the target `non-targeting`, and a targeting gRNA has another target.")))
    }
    list()
}

# Write a SingleCellExperiment as a data file.
#   assay "counts"          counts of the measured modality, features x cells
#   colData                 the cell covariates in one table, a modality's columns
#                           prefixed with its name (`rna_n_umis`, `grna_n_umis`)
#   rowData                 feature metadata
#   altExp "gRNA"           gRNA counts; its rowData is the gRNA table (grna_target, ...)
#   other altExps           further modalities over the same cells, written under
#                           their names
#   metadata(sce)           the entries of S10 that are known; `dataset_id` is kept
#                           if present; `main_modality` names the main
#                           experiment's modality ("rna" if absent); `title`
#                           and `source` (a URL, DOI or accession) are optional
fmt_write_screen <- function(sce, file, index_limit = NULL) {
    if (file.exists(file)) fmt_stop(file, " exists; a data file is written once.")
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    stamp <- S4Vectors::metadata(sce)
    main <- if (is.null(stamp$main_modality)) DEFAULT_MAIN_MODALITY else stamp$main_modality
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
    has_positions <- any(vapply(intersect(GRNA_COORDINATES, names(grna_table)), function(name) any(!is.na(grna_table[[name]])), logical(1)))
    if (has_positions && is.null(stamp$genome_build)) fmt_stop("The gRNA table has positions: metadata(sce)$genome_build must say which genome build they are in.")
    if (!is.null(grna_table$vector_id)) {
        if (anyNA(grna_table$vector_id)) fmt_stop("Every gRNA needs a `vector_id` when the column is present.")
        for (fault in vector_faults(rownames(grna), as.character(grna_table$vector_id), as.character(grna_table$grna_target))) fmt_stop(fault$message)
    }
    if (!is.null(grna_table$grna_type)) {
        for (fault in type_faults(as.character(grna_table$grna_type), as.character(grna_table$grna_target))) fmt_stop(fault$message)
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
            # scalar entries at the top level of `uns` (S10)
            metadata = c(
                list(
                    perturb_seq_format = scalar(FORMAT_VERSION),
                    dataset_id = scalar(if (is.null(stamp$dataset_id)) new_id() else stamp$dataset_id)
                ),
                # optional: a title for people, and what is known about the screen (S10.5)
                lapply(stamp[intersect(OPTIONAL_METADATA, names(stamp)[!vapply(stamp, is.null, logical(1))])], scalar)
            )
        ),
        file
    )
    for (name in names(experiments)) write_counts(file, name, experiments[[name]]$counts, index_limit)
    written <- TRUE
    invisible(file)
}

check_format <- function(format, what, known = FORMAT_VERSION) {
    if (is.null(format_parts(format))) fmt_stop(what, " has no format version.")
    if (!known_major(format, known)) {
        fmt_stop(what, " is in format ", format, "; this software knows format ", known, " and its minor versions.")
    }
}

# What software must refuse before using a data file it has read with
# MuData::readH5MU: a major format version it does not know (S10.3), and a gRNA
# table whose vectors break S8.7 or S8.8 (S8.9). Cheap: it looks at the metadata
# and the gRNA table only.
check_screen <- function(mae) {
    check_format(as.vector(S4Vectors::metadata(mae)$perturb_seq_format)[1], "The data file")
    grna_table <- SummarizedExperiment::rowData(mae[[GRNA]])
    if (!is.null(grna_table$vector_id)) {
        faults <- vector_faults(rownames(grna_table), as.character(grna_table$vector_id), as.character(grna_table$grna_target))
        for (fault in faults) fmt_stop("The data file: ", fault$message)
    }
}

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

# validate_screen() checks a data file against the numbered rules of SPEC.md.
# It reads the file directly, builds no object, and reports every problem it
# finds: a data frame with one row per problem (`rule`, `where`, `message`),
# with no rows when the file conforms. The Python library has the same function
# and reports the same rules.

problem_log <- function() {
    found <- list()
    list(
        add = function(rule, where, ...) {
            found[[length(found) + 1L]] <<- data.frame(rule = rule, where = where, message = paste0(...))
        },
        all = function() {
            if (length(found) == 0L) return(data.frame(rule = character(), where = character(), message = character()))
            out <- unique(do.call(rbind, found))
            out <- out[order(out$rule, out$where, out$message, method = "radix"), ]
            rownames(out) <- NULL
            out
        }
    )
}

MODALITY_NAME <- SCREEN_SPEC$modalities$name_pattern

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
        log$add("S10.1", where, "needs the scalar text entries `perturb_seq_format` and `dataset_id`.")
    } else if (is.null(format_parts(read(paste0(where, "/perturb_seq_format"))))) {
        log$add("S10.1", paste0(where, "/perturb_seq_format"), "is not a version of the form MAJOR.MINOR.")
    } else {
        if (!known_major(read(paste0(where, "/perturb_seq_format")))) log$add("S10.3", paste0(where, "/perturb_seq_format"), "is a format this validator does not know.")
        if (!nzchar(read(paste0(where, "/dataset_id")))) log$add("S10.1", paste0(where, "/dataset_id"), "is empty.")
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
                # not conforming; a reader still accepts any whole-number type here (S4.4)
                log$add(rules[[array]], paste0(x, "/", array), "is not stored as ", sub("int", "", kind), "-bit signed integers", why, ".")
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
                log$add("S5.3", path, "is not a valid column name (letters, digits, underscores and dots, starting with a letter).")
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
            if (any(given) && nrow(record) > 0L && nrow(entry("genome_build")) == 0L) log$add("S8.14", table, "The gRNA table has positions and the file has no `genome_build`.")
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
        # the kind of each gRNA, where recorded (S8.15, S8.16)
        if ("grna_type" %in% columns$name && "grna_target" %in% columns$name) {
            types <- tryCatch(text("grna_type"), error = function(e) NULL)
            targets <- tryCatch(text("grna_target"), error = function(e) NULL)
            if (!is.null(types) && !is.null(targets) && !anyNA(targets)) {
                for (fault in type_faults(types, targets)) log$add(fault$rule, paste0(table, "/grna_type"), fault$message)
            }
        }
    }
    log$all()
}
