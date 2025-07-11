#' The function spARI is the main function to calculate spRI and spARI values proposed
#'     in the paper.
#'
#' @param r_labels Annotated labels of all spots/cells. Can be numeric vector or character vector.
#' @param c_labels Estimated labels obtained by a certain spatial clustering method.
#' @param coords Spatial coordinates (2 columns). 1st column: first dimension coordinate.
#'     2nd column: second dimension coordinate. Default is NULL.
#' @param dist_mat Distance matrix provided by users. If both coords and dist_mat
#'     are provided, we will directly use the distance matrix. Default is NULL.
#' @param Is_spmat Logical; if TRUE, dist_mat is a sparse matrix ("dgCMatrix" or "dgTMatrix"). Default is FALSE.
#' @param alpha_val Coefficient belongs to the open interval (0, 1) to keep a positive gap
#'     between the maximal weight of the disagreement pair and the weight one of the agreement pair.
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
#' res_value1 = spARI(r_labels=true_labels, c_labels=c1_labels, coords)
#' res_value2 = spARI(r_labels=true_labels, c_labels=c2_labels, coords)
#' cat(paste0("1st method: spRI=", round(res_value1[1], 3), ", spARI=", round(res_value1[2], 3)))
#' cat(paste0("2nd method: spRI=", round(res_value2[1], 3), ", spARI=", round(res_value2[2], 3)))
#'
#' @export

spARI = function(r_labels, c_labels, coords=NULL, dist_mat=NULL, Is_spmat=FALSE, alpha_val=0.8) {

  if (length(unique(r_labels)) == 1 & length(unique(c_labels)) == 1) {
    return(c(1, 1))
  }

  n_objects = length(r_labels)
  n_obj_choose = choose(n_objects, 2)

  if (is.null(coords) & is.null(dist_mat)) {
    stop("Please provide either the spatial coordinates or the distance matrix!")
  } else if (!is.null(coords) & is.null(dist_mat)) {
    ## Coordinates normalization
    coords[,1] = (coords[,1] - min(coords[,1])) / (max(coords[,1]) - min(coords[,1]))
    coords[,2] = (coords[,2] - min(coords[,2])) / (max(coords[,2]) - min(coords[,2]))
    dist_mat = stats::dist(coords)
  }


  ## f function
  f_func = function(t) {
    return(alpha_val * exp(- t^2))
  }

  ## h function
  h_func = function(t) {
    return(alpha_val * (1 - exp(- t^2)))
  }

  ## f function (sparse mat)
  f_func_sp = function(t_sp) {
    return(alpha_val * exp(- t_sp@x^2))
  }

  ## h function (sparse mat)
  h_func_sp = function(t_sp) {
    return(alpha_val * (1 - exp(- t_sp@x^2)))
  }

  ## sum of f function
  sum_f_func = function(obj_pairs, dist_mat) {
    dist_vec = dist_mat[(n_objects - 0.5) * obj_pairs[,1] - 0.5 * obj_pairs[,1]^2 + obj_pairs[,2] - n_objects]
    f_vals = f_func(dist_vec)
    return(sum(f_vals))
  }

  ## sum of h function
  sum_h_func = function(obj_pairs, dist_mat) {
    dist_vec = dist_mat[(n_objects - 0.5) * obj_pairs[,1] - 0.5 * obj_pairs[,1]^2 + obj_pairs[,2] - n_objects]
    h_vals = h_func(dist_vec)
    return(sum(h_vals))
  }

  ## sum of f function (sparse mat)
  sum_f_func_sp = function(obj_pairs, dist_mat_sp) {
    dist_vec = dist_mat_sp[obj_pairs]
    f_vals = f_func(dist_vec)
    return(sum(f_vals))
  }

  ## sum of h function (sparse mat)
  sum_h_func_sp = function(obj_pairs, dist_mat_sp) {
    dist_vec = dist_mat_sp[obj_pairs]
    h_vals = h_func(dist_vec)
    return(sum(h_vals))
  }


  ## Compute spRI
  sg_pairs <- generate_sg_pairs_int(c_labels, r_labels) + 1
  gs_pairs <- generate_gs_pairs_int(c_labels, r_labels) + 1


  if (Is_spmat) {
    sums_f_and_h = sum_f_func_sp(sg_pairs, dist_mat) +
      sum_h_func_sp(gs_pairs, dist_mat)
  } else {
    sums_f_and_h = sum_f_func(sg_pairs, dist_mat) +
      sum_h_func(gs_pairs, dist_mat)
  }

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
    f_vals_total = f_func_sp(dist_mat)
    h_vals_total = h_func_sp(dist_mat)
  } else {
    f_vals_total = f_func(dist_mat)
    h_vals_total = h_func(dist_mat)
  }


  E_spRI = p*q*n_obj_choose + (1-p)*(1-q)*n_obj_choose +
    (1-p)*q*sum(f_vals_total) + p*(1-q)*sum(h_vals_total)
  E_spRI = E_spRI / n_obj_choose

  spARI_value = (spRI_value - E_spRI) / (1 - E_spRI)


  ## Outputs
  outputs = c(spRI_value, spARI_value)
  names(outputs) = c("spRI", "spARI")
  return(outputs)
}

