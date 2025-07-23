Changes in version 0.99.0 (2025-07-15)
+ Submitted to Bioconductor



Changes in version 0.99.1 (2025-07-23)

+ Generate spARI_example_data.rda in the directory `data/` and finish man documenation on the data object

+ Remove `exportPattern("^[[:alpha:]]+")` and export explicitly

+ Remove `vignette/spARI.html` 

+ Expand to work with a SpatialExperiment object. Add an argument `spe` as input of the functions `spARI` and `perm_test`

+ Add the following codes
  
  ```R
  r_labels = match(r_labels, unique(r_labels)) 
  c_labels = match(c_labels, unique(c_labels))
  ```
  
  in the `spARI` function to make it compatible with the character labels
  
+ Add the following codes

  ```R
  if (length(unique(r_labels)) == length(r_labels) & length(unique(c_labels)) == length(c_labels)) {
    cat("The spARI value is always equal to one!\n")
    return(invisible(NULL))
  }
  ```

  in the `perm_test` function to additionally consider the case that each single object forms a cluster
  
+ Change the code `p_value = sum(spARI_sim_record >= spARI_obs) / replicate_times` to `p_value = sum(spARI_sim_record > spARI_obs) / replicate_times` to align with the context in the paper