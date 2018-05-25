#' Creates a proper SingleCellExperiment object from expression tables and metadata
#'
#' @param data.files A named list with file paths to tab-separated expression matrices
#' @param metadata.files A list with file paths to tab-separated metadata matrices
#' @param cell.identifier The column name of cell identifiers in metadata.files (e.g. Sample, Cell)
#' @param gene.identifier The column name of gene identifieres in data.files (e.g. gene_id)
#' @param annotate.genes Boolean, whether or not bioMart should be used to annotate genes
#' @param biomart.dataset The dataset that is used to fetch annotations (e.g. mmusculus_gene_ensembl)
#' @param biomart.filter The filter that is used to match annotations to genes (e.g. ensembl_gene_id, if row.names are Ensembl Gene IDs)
#' @param biomart.fields A character vector of field names to fetch
#'
#' @return A SingleCellExperiment object.
#'
#' @export
createSingleCellExperiment <- function(data.files, metadata.files = NULL, cell.identifier = "Sample", gene.identifier = "gene_id", annotate.genes = F, biomart.dataset = "mmusculus_gene_ensembl", biomart.filter = "ensembl_gene_id", biomart.fields = c("ensembl_gene_id","external_gene_name","gene_biotype","description")) {
  #
  # Read expression data
  #
  expression <- lapply(data.files, function(p) data.frame(data.table::fread(input = p, sep = "\t", header = T), row.names = 1))

  #
  # Annotate gene names
  #
  if (annotate.genes) {
    identifier <- rownames(expression[[1]])
    bm <- biomaRt::getBM(biomart.fields, filters = biomart.filter, values = identifier, mart = biomaRt::useMart("ensembl", dataset = biomart.dataset))
    m <- match(identifier, bm[, biomart.filter])

    data <- data.frame(bm[stats::na.omit(m), ])
    rownames(data) <- data[, biomart.filter]
    colnames(data) <- biomart.fields

    rowData <- data[,c(2:ncol(data))]

    if (nrow(rowData) != nrow(expression[[1]])) {
      warning(paste("Could not annotate all", gene.identifier, "using biomaRt. Omitting unannotated genes."))
      expression <- lapply(expression, function(l) l[which(!is.na(m)), ])
    }
  } else {
    rowData <- NULL
  }

  #
  # Assemble metadata
  #
  if (!is.null(metadata.files)) {
    metadata_ <- lapply(metadata.files, function(p) metadata <- data.table::fread(input = p, sep = "\t", header = T, stringsAsFactors = T, na.strings = "NA"))
    metadata <- Reduce(function(a, b) dplyr::left_join(a, b, by = cell.identifier), metadata_)

    m <- match(colnames(expression[[1]]), make.names(metadata[, get("cell.identifier")]))
    colData <- as.data.frame(metadata[stats::na.omit(m), ])
    rownames(colData) <- make.names(colData[, cell.identifier])
  } else {
    colData <- NULL
  }

  # Convert to matrix
  expression <- lapply(expression, as.matrix)

  if (!is.null(colData) & !is.null(rowData)) {
    return(SingleCellExperiment::SingleCellExperiment(assays = expression, colData = colData, rowData = rowData))
  }
  if (is.null(colData) & is.null(rowData)) {
    return(SingleCellExperiment::SingleCellExperiment(assays = expression))
  }
  if (is.null(colData)) {
    return(SingleCellExperiment::SingleCellExperiment(assays = expression, rowData = rowData))
  }
  if (is.null(rowData)) {
    return(SingleCellExperiment::SingleCellExperiment(assays = expression, colData = colData))
  }
}
