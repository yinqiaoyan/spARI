% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spARI.R
\name{spARI}
\alias{spARI}
\title{The function spARI is the main function to calculate spRI and spARI values proposed
    in the paper.}
\usage{
spARI(
  r_labels,
  c_labels,
  coords = NULL,
  dist_mat = NULL,
  Is_spmat = FALSE,
  alpha_val = 0.8
)
}
\arguments{
\item{r_labels}{Annotated labels of all spots/cells. Can be numeric vector or character vector.}

\item{c_labels}{Estimated labels obtained by a certain spatial clustering method.}

\item{coords}{Spatial coordinates (2 columns). 1st column: first dimension coordinate.
2nd column: second dimension coordinate. Default is NULL.}

\item{dist_mat}{Distance matrix provided by users. If both coords and dist_mat
are provided, we will directly use the distance matrix. Default is NULL.}

\item{Is_spmat}{Logical; if TRUE, dist_mat is a sparse matrix ("dgCMatrix" or "dgTMatrix"). Default is FALSE.}

\item{alpha_val}{Coefficient belongs to the open interval (0, 1) to keep a positive gap
between the maximal weight of the disagreement pair and the weight one of the agreement pair.
Default is 0.8.}
}
\value{
spARI returns an R numeric including the following information.
\item{spRI_value}{numeric, the spRI value calculated by r_labels and c_labels}
\item{spARI_value}{numeric, the spARI value calculated by r_labels and c_labels}
}
\description{
The function spARI is the main function to calculate spRI and spARI values proposed
    in the paper.
}
\examples{
library(spARI)
### --- Import example data --- ###
# (1) true_labels: ground truth of 160 spots (input as r_labels)
# (2) c1_labels: one partition results of these spots (input as c_labels)
# (3) c2_labels: another partition results of these spots (input as c_labels)
# (4) coords: spatial coordinates of these spots (2 columns)
data(spARI_example_data)
### --- Compute spRI and spARI --- ###
res_value1 = spARI(r_labels=true_labels, c_labels=c1_labels, coords)
res_value2 = spARI(r_labels=true_labels, c_labels=c2_labels, coords)
cat(paste0("1st method: spRI=", round(res_value1[1], 3), ", spARI=", round(res_value1[2], 3)))
cat(paste0("2nd method: spRI=", round(res_value2[1], 3), ", spARI=", round(res_value2[2], 3)))

}
