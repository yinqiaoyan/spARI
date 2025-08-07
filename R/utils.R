#' @useDynLib spARI, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
NULL


#' Example data for spARI
#'
#' This dataset includes simulated spatial clustering results for demonstration purposes.
#' It contains the following objects:
#' \itemize{
#'   \item \code{true_labels}: Ground-truth labels of 160 spots.
#'   \item \code{c1_labels}: First clustering result to compare.
#'   \item \code{c2_labels}: Second clustering result to compare.
#'   \item \code{coords}: Spatial coordinates (matrix with 2 columns).
#' }
#'
#' @docType data
#' @keywords datasets
#' @name spARI_example_data
#' @usage data(spARI_example_data)
#' @format An environment containing 4 objects:
#' \describe{
#'   \item{true_labels}{Numeric or factor vector of length 160.}
#'   \item{c1_labels}{Numeric or factor vector of length 160.}
#'   \item{c2_labels}{Numeric or factor vector of length 160.}
#'   \item{coords}{Numeric matrix of size 160x2.}
#' }
NULL


#' Generate object pairs (sg pairs)
#'
#' This function identifies all object pairs that belong to the same cluster in the clustering partition
#' but different groups in the reference partition
#'
#' @name generate_sg_pairs_int
#'
#' @param c_labels An integer vector of clustering labels.
#' @param r_labels An integer vector of reference labels.
#'
#' @return An integer matrix with two columns, where each row is a pair of object indices.
#'
#' @examples
#' c_labels <- c(1,1,2,2,2,3,3)
#' r_labels <- c(1,1,1,2,2,3,1)
#' generate_sg_pairs_int(c_labels, r_labels)
#'
#' @export
NULL


#' Generate object pairs (gs pairs)
#'
#' This function identifies all object pairs that belong to the same group in the reference partition
#' but different clusters in the clustering partition
#'
#' @name generate_gs_pairs_int
#'
#' @param c_labels An integer vector of clustering labels.
#' @param r_labels An integer vector of reference labels.
#'
#' @return An integer matrix with two columns, where each row is a pair of object indices.
#'
#' @examples
#' c_labels <- c(1,1,2,2,2,3,3)
#' r_labels <- c(1,1,1,2,2,3,1)
#' generate_gs_pairs_int(c_labels, r_labels)
#'
#' @export
NULL




