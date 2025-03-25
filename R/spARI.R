#' The function spARI is the main function to calculate the spARI value proposed
#'     in this paper, which is the adjusted version of spRI.
#'
#' @param r_labels Annotated labels of all spots/cells. Can be numeric vector or character vector.
#' @param c_labels Estimated labels obtained by a certain spatial clustering method.
#' @param coords Spatial coordinates (2 columns). 1st column: first dimension coordinate.
#'     2nd column: second dimension coordinate.
#' @param alpha_val Coefficient belongs to the open interval (0, 1) to keep a positive gap
#'     between the maximal weight of the disagreement pair and the weight one of the agreement pair.
#'     Default is 0.8.
#' @param print_time Logical; if TRUE, the total execution time is printed. Default is FALSE.
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
#' res_value1 = spARI(r_labels=true_labels, c_labels=c1_labels, coords, print_time = TRUE)
#' res_value2 = spARI(r_labels=true_labels, c_labels=c2_labels, coords, print_time = TRUE)
#' cat(paste0("1st method: spRI=", round(res_value1[1], 3), ", spARI=", round(res_value1[2], 3)))
#' cat(paste0("2nd method: spRI=", round(res_value2[1], 3), ", spARI=", round(res_value2[2], 3)))
#'
#' @export

spARI = function(r_labels, c_labels, coords, alpha_val=0.8, print_time=FALSE) {

  if (length(unique(r_labels)) == 1 & length(unique(c_labels)) == 1) {
    return(c(1, 1))
  }

  n_objects = length(r_labels)
  n_obj_choose = choose(n_objects, 2)

  ## Coordinates normalization
  coords[,1] = (coords[,1] - min(coords[,1])) / (max(coords[,1]) - min(coords[,1]))
  coords[,2] = (coords[,2] - min(coords[,2])) / (max(coords[,2]) - min(coords[,2]))

  dist_mat = stats::dist(coords)

  ## f function
  f_func = function(t) {
    return(alpha_val * exp(- t^2))
  }

  ## h function
  h_func = function(t) {
    return(alpha_val * (1 - exp(- t^2)))
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


  ## Compute spRI
  cat("=== Computing spRI ===\n")
  s_time = Sys.time()
  same_group_c_pairs <- outer(c_labels, c_labels, `==`)
  same_group_r_pairs <- outer(r_labels, r_labels, `==`)
  sg_pairs = which(same_group_c_pairs & !same_group_r_pairs, arr.ind = TRUE)  # s*g type
  gs_pairs = which(same_group_r_pairs & !same_group_c_pairs, arr.ind = TRUE)  # g*s type

  sg_pairs = sg_pairs[sg_pairs[,1] < sg_pairs[,2], ]
  gs_pairs = gs_pairs[gs_pairs[,1] < gs_pairs[,2], ]

  # dist_vec_sg = dist_mat[(n_objects - 0.5) * sg_pairs[,1] - 0.5 * sg_pairs[,1]^2 + sg_pairs[,2] - n_objects]
  # dist_vec_gs = dist_mat[(n_objects - 0.5) * gs_pairs[,1] - 0.5 * gs_pairs[,1]^2 + gs_pairs[,2] - n_objects]

  sums_f_and_h = sum_f_func(sg_pairs, dist_mat) +
    sum_h_func(gs_pairs, dist_mat)
  # sums_f_and_h = sum(f_func(dist_vec_sg)) +
  #   sum(h_func(dist_vec_gs))


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
  cat("=== Computing spARI ===\n")
  p = 0.5 * (parSum_r_sqsum - n_objects) / n_obj_choose
  q = 0.5 * (parSum_c_sqsum - n_objects) / n_obj_choose

  f_vals_total = f_func(dist_mat)
  h_vals_total = h_func(dist_mat)

  E_spRI = p*q*n_obj_choose + (1-p)*(1-q)*n_obj_choose +
    (1-p)*q*sum(f_vals_total) + p*(1-q)*sum(h_vals_total)
  # E_spRI = (p*q + (1-p)*(1-q))*n_obj_choose +
  #   (1-p)*q*sum(f_vals_total) + p*(1-q)*sum(h_vals_total)  # slower !?
  E_spRI = E_spRI / n_obj_choose

  spARI_value = (spRI_value - E_spRI) / (1 - E_spRI)

  e_time = Sys.time()
  exeTime = e_time - s_time

  if (print_time) {
    print(exeTime)
  }

  outputs = c(spRI_value, spARI_value)
  return(outputs)
}

