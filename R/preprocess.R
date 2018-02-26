#' Generate a mutationhash for a vector of barcodes.
#'
#' @param barcodes Character vector of barcodes
#' @param edit_distance Maximum allowed edit distance
#'
#' @return A data.frame with two columns: A column containing barcodes and a column containing possible substituted barcodes.
#'
#' @export
generate_mutationhash <- function(barcodes, edit_distance = 1) {
  muthash <- data.frame(barcode = barcodes, mutation = barcodes)
  for (barcode in barcodes) {
    set <- substitution_set(barcode, edit_distance)

    df <- data.frame(barcode = rep(barcode, length(set)), mutation = set)
    muthash <- rbind(muthash, df)
  }
  unique(muthash)
}

#' Generate an substitution set of all strings that are \code{nedit} edits away from \code{str} under a given alphabet.
#'
#' @param str The string that is substituted.
#' @param nedit Maximum allowed edit distance.
#' @param alphabet Possible alphabet, from which substitutions are drawn.
#'
#' @return A vector of strings.
#'
substitution_set <- function(str, nedit, alphabet = c("A", "T", "G", "C")) {
  strlen <- stringr::str_length(str)
  indices <- utils::combn(seq_len(strlen), nedit, simplify = F)
  substitutions = expand.grid(list(alphabet)[rep(1,nedit)], stringsAsFactors = F)

  strlist <- unlist(strsplit(str, fixed = T, split = ""))
  unique(unlist(apply(substitutions, 1, function(s) {
    set <- lapply(indices, function(j) {
      strlist[j] <- s
      return(stringr::str_c(unlist(strlist), collapse = ""))
    })
  })))
}

