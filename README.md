## spARI

Current version: 0.99.0 (20250715)

The R package **spARI** is designed to compute two novel clustering evaluation metrics—the **spatially aware Rand index (spRI)** and its adjusted version (**spARI**)—for assessing spatial transcriptomics clustering. Unlike the traditional Rand index (RI) and adjusted Rand index (ARI), spRI and spARI incorporate spatial distance information into clustering evaluation. When comparing two partitions, spRI assigns a weight to each disagreement pair—two objects in the same cluster of one partition but in different clusters of the other—based on their spatial distance, allowing for a more refined distinction between disagreement pairs. The spRI value ranges between zero and one, while the spARI value is less than one with an expected value of zero. Higher spARI values indicate greater clustering accuracy and a more compact spatial structure. 

To implement this package, users are required to input the two partitions, and to optionally provide either the spatial coordinates or the distance matrix. The main function "spARI" then automatically computes spRI and spARI in a fast and user-friendly manner. If spatial coordinates are provided, a built-in normalization step is applied to remove the effect of spatial location unit. For large-scale spatial data (e.g., 100,000 objects), users can input a sparse distance matrix of class "dgCMatrix" or "dgTMatrix", for instance, by retaining only distances to the $k$ nearest neighbors for each object. The `spARI` package is compatible with Windows, Linux, and macOS, and can be easily installed on all three platforms.



## Prerequisites and Installation

1. R version >= 4.1.3
2. CRAN package: Rcpp, stats, parallel, Matrix
3. Install the package `spARI`

```R
devtools::install_github("yinqiaoyan/spARI")
```



## Datasets information

The data description is given in the following table.

|  ST Dataset  | Spot/Cell number | Gene number |                        Download links                        |
| :----------: | :--------------: | :---------: | :----------------------------------------------------------: |
| DLPFC 151509 |      4,789       |   33,538    |        https://github.com/LieberInstitute/spatialLIBD        |
|   STARmap*   |      1,207       |    1,020    | Raw data: http://sdmbench.drai.cn/tcm/download/?file_path=/mnt/JINGD/data/file/sdmbench/db/STARmap_20180505_BY3_1k.h5ad  <br/>Cell type annotation: https://drive.google.com/drive/folders/1I1nxheWlc2RXSdiv24dex3YRaEh780my?usp=sharing |



## Example Code

### 1. Compute spRI and spARI using spatial coordinates

The following code shows an example (the third simulation study in the manuscript) that runs the main function "spARI" in our package.

Import the required R package.

```R
library(spARI)
```

Load the example data stored in this package. The example data includes:

* true_labels: the reference partition
* c1_labels: the first clustering partition
* c2_labels: the second clustering partition
* coords: spatial coordinates of all the objects

```R
data("spARI_example_data")
```

Run "spARI" function to compute spRI and spARI. The meaning of each argument in the main function is listed below.

* r_labels: the reference partition
* c_labels: the clustering partition
* coords: spatial coordinates of all the objects. Default is NULL
* dist_mat: distance matrix calculated by users. Default is NULL
* Is_spmat: a logical value indicating whether "dist_mat" is a sparse matrix. Default is FALSE
* alpha_val: coefficient belongs to the open interval (0, 1) to keep a positive gap between the maximal weight of the disagreement pair and the weight one of the agreement pair. Default is 0.8

```R
res_value1 = spARI(r_labels=true_labels, c_labels=c1_labels, coords=coords)
res_value2 = spARI(r_labels=true_labels, c_labels=c2_labels, coords=coords)
```

Users can simply run the code `example("spARI")` to carry out this example.

```R
library(BACT)
example("spARI")
# 1st clustering partition: spRI=0.926, spARI=0.672
# 2nd clustering partition: spRI=0.917, spARI=0.632
```

### 2. Compute spRI and spARI using sparse distance matrix

In this example, we generate two partitions for 100,000 objects with simulated spatial coordinates to examine the computational efficiency of the `spARI` package for large-scale spatial data.

Since storing a dense matrix of 100,000$\times$100,000 would require more than 16GB of memory, we generate a sparse distance matrix of the objects. 

```R
## Install FNN package
if (!require("FNN", quietly = TRUE))
  install.packages("FNN")

library(FNN)
library(Matrix)

## Define sparse distance matrix generation function
build_symmetric_knn_distance_matrix <- function(coord_mat, k = 5) {
  N <- nrow(coord_mat)
  
  # Find k nearest neighbors for each object (exclude itself)
  knn_res <- FNN::get.knn(coord_mat, k = k)
  i_idx <- rep(1:N, each = k)
  j_idx <- as.vector(t(knn_res$nn.index))
  d_val <- as.vector(t(knn_res$nn.dist))
  
  # Generate sparse matrix (symmetric)
  D_temp <- sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = d_val,
    dims = c(N, N)
  )
  
  D1 <- as(D_temp, "dgTMatrix")
  D2 <- as(t(D_temp), "dgTMatrix")
  
  i_all <- c(D1@i, D2@i)
  j_all <- c(D1@j, D2@j)
  x_all <- c(D1@x, D2@x)
  
  # Combine and aggregate max
  df <- data.frame(i = i_all + 1, j = j_all + 1, x = x_all)  # +1 for R indexing
  df_agg <- aggregate(x ~ i + j, data = df, FUN = max)
  
  # Create symmetric sparse matrix
  D_sym <- sparseMatrix(i = df_agg$j, j = df_agg$i, x = df_agg$x)
  
  return(D_sym)
}

## Generate sparse distance matrix
set.seed(123)
N <- 1e5  # 100,000 objects
coords <- matrix(runif(N*2), N, 2)  # coordinates ranging between 0 and 1
dist_mat <- build_symmetric_knn_distance_matrix(coords, k = 5)

dim(dist_mat)
# 100000 100000
length(dist_mat@x)
# 595850  # number of non-zero elements (0.006% of total N^2 elements)
```

**Remark.** For each object, we initially retain only the distances to its $k=5$ nearest neighbors (recorded entries in the distance matrix). However, since the $k$-nearest neighbor relationship is not symmetric (i.e., object $i$ may be among the $k$ nearest neighbors of object $j$, but $j$ is not necessarily among the $k$ nearest neighbors of $i$), the resulting matrix $D_{\text{temp}}$ is asymmetric. To ensure symmetry, we include the distance between two objects if either one is among the $k$ nearest neighbors of the other.

Then we generate the synthetic reference and clustering partitions, with the clustering result containing 10,000 misclustered objects.

```R
K = 15
set.seed(123)

# reference partition
true_labels = sample(1:K, N, replace = T)
# clustering partition
c_labels = true_labels
ids = sample(which(true_labels == 1), 5000)
c_labels[ids] = 2
ids = sample(which(true_labels == 3), 5000)
c_labels[ids] = 4
```

Compute spRI and spARI values.

```R
library(spARI)

stime = Sys.time()
res = spARI(true_labels, c_labels, dist_mat=dist_mat)
etime = Sys.time()
print(res)
#      spRI     spARI 
# 0.9940251 0.9550960 
print(etime-stime)
# Time difference of 5.503251 secs
```

### 3. Conducting hypothesis testing

We conduct the permutation test for the observed spARI value of the two clustering partitions in the example data. 

First, we import `spARI` package and load the example data

```R
library(spARI)
data("spARI_example_data")
```

Then we carry out the permutation test for the clustering partitions "c1_labels" and "c2_labels".

```R
perm_test(r_labels=true_labels, c_labels=c1_labels, coords=coords)
# spARI_obs   p-value 
# 0.6716089 0.0000000

perm_test(r_labels=true_labels, c_labels=c2_labels, coords=coords)
# spARI_obs   p-value 
#  0.631726  0.000000 
```

Both the p-values are smaller than a given significance level (e.g., 0.01), indicating that the observed spARI value is significantly larger than zero, and the null hypothesis "spARI=0" is rejected.



## Contact

If you have any questions regarding this package, please contact Yinqiao Yan at [yinqiaoyan@bjut.edu.cn](mailto:yinqiaoyan@bjut.edu.cn).

