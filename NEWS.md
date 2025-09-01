# spARI 0.99.12 (2025-09-01)

- Updated the README file.

# spARI 0.99.11 (2025-08-30)

- Removed `parallel` from the `Imports` field in the DESCRIPTION file.

# spARI 0.99.10 (2025-08-30)

- Deleted system file `.Rapp.history`

# spARI 0.99.9 (2025-08-30)

- Deleted system files

# spARI 0.99.8 (2025-08-30)

- Updated the README file
- Updated the descriptions of parameters and examples in R functions `spARI` and `perm_test` 
- Updated the code `sum_f_vals_total = sum(f_func(dist_mat_Lvec))` to `sum_f_vals_total = sum(f_func(dist_mat_Lvec[dist_mat_Lvec != 0]))`
- Added `set.seed()` to the example of `perm_test` when `use_parallel = FALSE`
- Replaced all occurrences of `=` used for assignment (except in function arguments) with `<-`

# spARI 0.99.7 (2025-08-07)

- Converted `spARI_example_data` into a list object

# spARI 0.99.6 (2025-08-07)

- Updated `utils.R` file

# spARI 0.99.5 (2025-08-07)

- Addressed code styling suggestions (replaced `1:n` by `seq_len()`, replaced `sapply` by `vapply`)
- Removed usage of `set.seed()` in package code
- Added missing `\value{}` sections in man pages

# spARI 0.99.4 (2025-08-07)

- Removed `.DS_Store` files in the package

# spARI 0.99.3 (2025-08-07)

- Updated the Readme file

# spARI 0.99.2 (2025-08-07)

- Subscribed to the Bioc-devel mailing list

# spARI 0.99.1 (2025-07-23)

- Generated spARI_example_data.rda in the directory `data/` and finished man documenation on the data object
- Removed `exportPattern("^[[:alpha:]]+")` and export explicitly
- Removed `vignette/spARI.html` 
- Expanded support for `SpatialExperiment` objects; added `spe` argument to `spARI` and `perm_test` functions
- Made `spARI` compatible with character labels by mapping `r_labels` and `c_labels` to numeric factors
- Added specific case handling in `perm_test` for singleton clusters
- Changed p-value calculation from `>=` to `>` in `perm_test`, aligning with the paper

# spARI 0.99.0 (2025-07-15)

- Submitted to Bioconductor