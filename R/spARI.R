#' The function spARI is the main function to calculate spRI and spARI values proposed
#'     in the paper.
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
#'
#' @return spARI returns an R numeric including the following information.
#' \item{spRI_value}{numeric, the spRI value calculated by r_labels and c_labels}
#' \item{spARI_value}{numeric, the spARI value calculated by r_labels and c_labels}
#'
#' @examples
#' library(spARI)
#' ### --- Import example data --- ###
#' # (1) true_labels: ground truth of 160 spots (input as r_labels)
#' # (2) c1_labels: one partition results of these spots (input as c_labels)
#' # (3) c2_labels: another partition results of these spots (input as c_labels)
#' # (4) coords: spatial coordinates of these spots (2 columns)
#' data(spARI_example_data)
#' ### --- Compute spRI and spARI --- ###
#' res1 = spARI(r_labels=true_labels, c_labels=c1_labels, coords=coords)
#' res2 = spARI(r_labels=true_labels, c_labels=c2_labels, coords=coords)
#' cat(paste0("1st: spRI=", round(res1[1], 3), ", spARI=", round(res1[2], 3), "\n"))
#' cat(paste0("2nd: spRI=", round(res2[1], 3), ", spARI=", round(res2[2], 3), "\n"))
#'
#' @export

spARI = function(r_labels, c_labels, coords=NULL, dist_mat=NULL,
                 f_func_input=NULL, h_func_input=NULL, alpha_val=0.8) {

  if (length(r_labels) == length(c_labels) &
      length(unique(r_labels)) == 1 &
      length(unique(c_labels)) == 1) {
    ## Outputs
    outputs = c(1, 1)
    names(outputs) = c("spRI", "spARI")
    return(outputs)
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
    Is_spmat = FALSE
  } else if (length(r_labels) != length(c_labels) |
             length(r_labels) != nrow(dist_mat) |
             length(c_labels) != nrow(dist_mat)) {
    stop("Please align the length of r_label, the length of c_label, and the dimension of dist_mat!")
  } else {
    Is_spmat = ("dgCMatrix" %in% class(dist_mat) | "dgTMatrix" %in% class(dist_mat))
  }

  n_objects = length(r_labels)
  n_obj_choose = choose(n_objects, 2)


  ## f function
  if (is.null(f_func_input)) {
    f_func = function(t) {
      return(alpha_val * exp(- t^2))
    }
  } else {
    f_func = f_func_input
  }

  ## h function
  if (is.null(h_func_input)) {
    h_func = function(t) {
      return(alpha_val * (1 - exp(- t^2)))
    }
  } else {
    h_func = h_func_input
  }


  ## Compute spRI
  sg_pairs <- generate_sg_pairs_int(c_labels, r_labels) + 1
  gs_pairs <- generate_gs_pairs_int(c_labels, r_labels) + 1

  dist_vec_sg = dist_mat[sg_pairs]
  dist_vec_gs = dist_mat[gs_pairs]

  sums_f_and_h = sum(f_func(dist_vec_sg[dist_vec_sg != 0])) +
    sum(h_func(dist_vec_gs[dist_vec_gs != 0]))

  unique_value_r = sort(unique(r_labels))
  unique_value_c = sort(unique(c_labels))
  K_R = length(unique_value_r)
  K_C = length(unique_value_c)

  nij_record = matrix(NA, nrow = K_R, ncol = K_C)
  for (rr in 1:K_R) {
    for (cc in 1:K_C) {
      tmp_rr = r_labels == unique_value_r[rr]
      tmp_cc = c_labels == unique_value_c[cc]
      nij_record[rr, cc] = sum(tmp_rr & tmp_cc)
    }
  }

  parSum_r = rowSums(nij_record)
  parSum_c = colSums(nij_record)
  parSum_r_sqsum = sum(parSum_r^2)
  parSum_c_sqsum = sum(parSum_c^2)

  num_A = n_obj_choose + sum(as.numeric(nij_record)^2) - 0.5 * (parSum_r_sqsum + parSum_c_sqsum)
  spRI_value = (num_A + sums_f_and_h) / n_obj_choose


  ## Compute spARI
  p = 0.5 * (parSum_r_sqsum - n_objects) / n_obj_choose
  q = 0.5 * (parSum_c_sqsum - n_objects) / n_obj_choose

  if (Is_spmat) {
    nonzero_vec = Matrix::tril(dist_mat)@x
    sum_f_vals_total = sum(f_func(nonzero_vec))
    sum_h_vals_total = sum(h_func(nonzero_vec))
  } else {
    dist_mat_Lvec = dist_mat[lower.tri(dist_mat, diag = FALSE)]
    sum_f_vals_total = sum(f_func(dist_mat_Lvec))
    sum_h_vals_total = sum(h_func(dist_mat_Lvec))
  }


  E_spRI = p*q + (1-p)*(1-q) +
    (1-p)*q*sum_f_vals_total / n_obj_choose +
    p*(1-q)*sum_h_vals_total / n_obj_choose

  if (E_spRI == 1) {
    ## Outputs
    outputs = c(1, 1)
    names(outputs) = c("spRI", "spARI")
    return(outputs)
  }

  spARI_value = (spRI_value - E_spRI) / (1 - E_spRI)


  ## Outputs
  outputs = c(spRI_value, spARI_value)
  names(outputs) = c("spRI", "spARI")
  return(outputs)
}

