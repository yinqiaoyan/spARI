test_that("spARI runs correctly and returns expected structure", {
  # load data
  data("spARI_example_data", package = "spARI")
  coords_norm = coords
  coords_norm[,1] = (coords[,1] - min(coords[,1])) / (max(coords[,1]) - min(coords[,1]))
  coords_norm[,2] = (coords[,2] - min(coords[,2])) / (max(coords[,2]) - min(coords[,2]))
  dist_mat = as.matrix(stats::dist(coords_norm))

  # basic call with coordinates
  result_coord <- spARI(true_labels, c1_labels, coords=coords)

  # basic call with distance matrix
  result_dist <- spARI(true_labels, c1_labels, dist_mat=dist_mat)

  # permutation test
  test_result <- perm_test(true_labels, c1_labels, coords=coords, use_parallel=FALSE)

  # check length of result
  expect_length(result_coord, 2)
  expect_length(result_dist, 2)
  expect_length(test_result, 2)

  # check names
  expect_named(result_coord, c("spRI", "spARI"))
  expect_named(result_dist, c("spRI", "spARI"))
  expect_named(test_result, c("spARI_obs", "p-value"))

  # check numeric
  expect_type(result_coord[["spRI"]], "double")
  expect_type(result_coord[["spARI"]], "double")

  expect_type(result_dist[["spRI"]], "double")
  expect_type(result_dist[["spARI"]], "double")

  expect_type(test_result[["spARI_obs"]], "double")
  expect_type(test_result[["p-value"]], "double")

  # check values are in [0, 1]
  expect_gte(result_coord[["spRI"]], 0)
  expect_lte(result_coord[["spRI"]], 1)
  expect_gte(result_coord[["spARI"]], 0)
  expect_lte(result_coord[["spARI"]], 1)

  expect_gte(result_dist[["spRI"]], 0)
  expect_lte(result_dist[["spRI"]], 1)
  expect_gte(result_dist[["spARI"]], 0)
  expect_lte(result_dist[["spARI"]], 1)

  expect_gte(test_result[["spARI_obs"]], 0)
  expect_lte(test_result[["spARI_obs"]], 1)
  expect_gte(test_result[["p-value"]], 0)
  expect_lte(test_result[["p-value"]], 1)
})
