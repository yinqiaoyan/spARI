#' The function perm_test is the main function to carry out the hypothesis test procedure
#'     proposed in the paper.
#'
#' @param r_labels Annotated labels of all spots/cells. Can be numeric vector or character vector.
#' @param c_labels Estimated labels obtained by a certain spatial clustering method.
#' @param coords Spatial coordinates (2 columns). 1st column: first dimension coordinate.
#'     2nd column: second dimension coordinate. Default is NULL.
#' @param dist_mat Distance matrix provided by users. If both coords and dist_mat
#'     are provided, we will directly use the distance matrix. Default is NULL.
#'     Please notice that if dist_mat is sparse, the weight function for object pairs without recorded distances
#'     degenerates to the setting used in the classical Rand index.
#' @param f_func_input R function; function f provided by users.
#' @param h_func_input R function; function h provided by users.
#' @param alpha_val Parameter in the default functions f and h, which belongs to the open interval (0, 1)
#'     to keep a positive gap between the maximal weight of the disagreement pair and the weight one of the agreement pair.
#'     Default is 0.8.
#' @param use_parallel Logical; if TRUE, use parallel code to permute the two partitions. Default is TRUE.
#' @param replicate_times Number of permutations for both the reference and clustering partitions. Default is 100.
#' @param random_seed Random seed for reproducibility. Default is 42.
#' @param spe SpatialExperiment object; stores various components of spatial transcriptomics data, including
#' spatialCoords: A matrix containing the spatial coordinates;
#' colData$cell_type: Annotated cell type labels for each spot or cell;
#' colData$cluster: Clustering labels for each spot or cell.
#' Default is NULL.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
#'
#' @return spARI returns an R numeric including the following information.
#' \item{spARI_obs}{numeric, the observed spARI value calculated by r_labels and c_labels}
#' \item{p_value}{numeric, the p-value of the hypothesis testing}
#'
#' @examples
#' library(spARI)
#' data(spARI_example_data)
#' true_labels = spARI_example_data$true_labels
#' c1_labels = spARI_example_data$c1_labels
#' c2_labels = spARI_example_data$c2_labels
#' coords = spARI_example_data$coords
#' test_res1 = perm_test(r_labels=true_labels, c_labels=c1_labels, coords=coords,
#'                       use_parallel=FALSE)
#' test_res2 = perm_test(r_labels=true_labels, c_labels=c2_labels, coords=coords,
#'                       use_parallel=FALSE)
#'
#' @export

perm_test = function(r_labels, c_labels, coords=NULL, dist_mat=NULL,
                     f_func_input=NULL, h_func_input=NULL, alpha_val=0.8,
                     use_parallel=TRUE, replicate_times=100, random_seed=42,
                     spe=NULL) {

  if (missing(r_labels) & missing(c_labels) & !is.null(spe)) {
    ## ST data is input in SpatialExperiment format
    if (is.null(colData(spe)$cell_type)) {
      stop("Please store cell type labels in \"colData$cell_type\"")
    }
    if (is.null(colData(spe)$cluster)) {
      stop("Please store clustering labels in \"colData$cluster\"")
    }
    coords = spatialCoords(spe)
    r_labels = colData(spe)$cell_type
    c_labels = colData(spe)$cluster
  }

  if (is.null(coords) & is.null(dist_mat)) {
    stop("Please provide either the spatial coordinates or the distance matrix!")
  } else if (!is.null(coords) & is.null(dist_mat)) {
    if (length(r_labels) != length(c_labels) |
        length(r_labels) != nrow(coords) |
        length(c_labels) != nrow(coords)) {
      stop("Please align the length of r_label, the length of c_label, and the dimension of coords!")
    }
    ## Coordinates normalization
    coords[,1] = (coords[,1] - min(coords[,1])) / (max(coords[,1]) - min(coords[,1]))
    coords[,2] = (coords[,2] - min(coords[,2])) / (max(coords[,2]) - min(coords[,2]))
    dist_mat = as.matrix(stats::dist(coords))
  } else if (length(r_labels) != length(c_labels) |
             length(r_labels) != nrow(dist_mat) |
             length(c_labels) != nrow(dist_mat)) {
    stop("Please align the length of r_label, the length of c_label, and the dimension of dist_mat!")
  }

  if ((length(unique(r_labels)) == length(r_labels) & length(unique(c_labels)) == 1) |
      (length(unique(c_labels)) == length(c_labels) & length(unique(r_labels)) == 1)) {
    cat("The spARI value is always equal to zero!\n")
    return(invisible(NULL))
  }


  # observed spARI value
  spARI_obs = spARI(r_labels, c_labels, dist_mat=dist_mat,
                    f_func_input = f_func_input, h_func_input = h_func_input,
                    alpha_val = alpha_val)[2]

  # permutation test (in parallel)
  if ((length(unique(r_labels)) == 1 & length(unique(c_labels)) == 1) |
      (length(unique(r_labels)) == length(r_labels) & length(unique(c_labels)) == length(c_labels))) {
    cat("The spARI value is always equal to one!\n")
    return(invisible(NULL))
  } else if (use_parallel) {
    ncores <- min(6, max(1, parallel::detectCores() - 2))
    cl <- parallel::makeCluster(ncores)
    parallel::clusterSetRNGStream(cl, iseed = random_seed)

    spARI_sim_record <- parallel::parLapply(
      cl, seq_len(replicate_times),
      fun = function(i, r_labels, c_labels, dist_mat, f_func_input, h_func_input, alpha_val) {
        r_labels_perm <- sample(r_labels)
        c_labels_perm <- sample(c_labels)
        spARI(r_labels_perm, c_labels_perm, dist_mat = dist_mat,
              f_func_input = f_func_input, h_func_input = h_func_input,
              alpha_val = alpha_val)[2]
      },
      r_labels, c_labels, dist_mat, f_func_input, h_func_input, alpha_val
    )

    parallel::stopCluster(cl)
    spARI_sim_record <- unlist(spARI_sim_record)
  } else {
    # set.seed(random_seed)
    spARI_sim_record <- vapply(seq_len(replicate_times), function(i) {
      r_labels_perm <- sample(r_labels)
      c_labels_perm <- sample(c_labels)
      spARI_val = spARI(r_labels_perm, c_labels_perm, dist_mat = dist_mat,
                        f_func_input = f_func_input, h_func_input = h_func_input,
                        alpha_val = alpha_val)[2]
      return(spARI_val)
    }, FUN.VALUE = numeric(1))
  }

  # compute p value
  p_value = sum(spARI_sim_record > spARI_obs) / replicate_times

  ## Outputs
  outputs = c(spARI_obs, p_value)
  names(outputs) = c("spARI_obs", "p-value")
  return(outputs)
}


