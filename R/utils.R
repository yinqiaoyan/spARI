#' @useDynLib spARI, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#'
NULL


#' Data: spARI_example_data
#'
#' spARI_example_data contains coords and gene_data_pc
#' true_labels: ground truth of 160 spots (input as r_labels)
#' c1_labels: one partition results of these spots (input as c_labels)
#' c2_labels: another partition results of these spots (input as c_labels)
#' coords: spatial coordinates of these spots (2 columns)
#' @name spARI_example_data
#' @docType data
#' @usage data(spARI_example_data)
NULL



