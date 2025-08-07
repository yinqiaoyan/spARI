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