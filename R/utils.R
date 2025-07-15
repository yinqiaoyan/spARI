#' @useDynLib spARI, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#'
NULL


#' Data: spARI_example_data
#'
#' spARI_example_data contains the following items:
#' @format A RData file:
#' \describe{
#'   \item{true_labels}{ground truth of 160 spots (input as r_labels)}
#'   \item{c1_labels}{one partition results of these spots (input as c_labels)}
#'   \item{c2_labels}{another partition results of these spots (input as c_labels)}
#'   \item{coords}{spatial coordinates of these spots (2 columns)}
#' }
#'
#' @name spARI_example_data
#' @docType data
#' @usage data(spARI_example_data)
NULL


#' Coordinates of Simulated Objects
#'
#' A simulated dataset containing spatial coordinates.
#'
#' @format A data frame with 2 columns:
#' \describe{
#'   \item{X_coord}{x-coordinate}
#'   \item{Y_coord}{y-coordinate}
#' }
#' @name coords
NULL

#' True Labels of Simulated Objects
#'
#' A manually generated true labels.
#' @format An integer vector.
#'
#' @name true_labels
NULL


#' Clustering Labels of Simulated Objects
#'
#' A manually generated clustering labels.
#' @format An integer vector.
#'
#' @name c1_labels
NULL


#' Clustering Labels of Simulated Objects
#'
#' Another manually generated clustering labels.
#' @format An integer vector.
#'
#' @name c2_labels
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




