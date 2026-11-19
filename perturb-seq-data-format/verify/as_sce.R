# For the checks only. MuData::readH5MU returns a MultiAssayExperiment and the
# prototype's R writer takes a SingleCellExperiment; this converts the first to
# the second. `response_modality` picks the modality that becomes the main
# experiment; it may be left out when the file has one modality besides `grna`.
#     source("verify/as_sce.R")   # after format/setup.R
read_as_sce <- function(file, response_modality = NULL) {
    mae <- withCallingHandlers(
        MuData::readH5MU(file),
        warning = function(w) {
            if (grepl("sampleMap", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
        }
    )
    stamp <- S4Vectors::metadata(mae)
    stamp <- lapply(stamp[intersect(names(SCREEN_SPEC$metadata$entries), names(stamp))], function(x) as.vector(x)[1])
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
    S4Vectors::metadata(sce) <- c(
        list(perturb_seq_format = stamp$perturb_seq_format, dataset_id = stamp$dataset_id, response_modality = main),
        stamp[intersect(OPTIONAL_METADATA, names(stamp))]
    )
    sce
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
