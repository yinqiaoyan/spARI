## spARI

Current `spARI` R package: version 1.0 

The R package `spARI` is built to the spatially aware Rand index (spRI) as well as its adjusted version (spARI) that incorporate the spatial distance information into the clustering evaluation. When comparing two partitions, spRI provides a pair of objects that are in the same cluster of one partition but are in different clusters of the other partition (called disagreement pair) with a weight relying on the distance of the two objects, and spARI is an adjustment of spRI that corrects for random chances such that its expectation takes on the zero value under an appropriate random null model. This R package can be installed in Windows, Linux, and Mac OS.



## Prerequisites and Installation

1. R version >= 4.1.3.
2. CRAN package: stats (>=4.1.3)
3. Install the package `spARI`.

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

The following code shows an example (the third simulation study in the manuscript) that runs the main function "spARI" in our package.

Import the required R packages.

```R
library(spARI)
```

 Read the example data stored in this package. The example data includes:

* true_labels: the reference partition
* c1_labels: the first clustering partition
* c2_labels: the second clustering partition
* coords: spatial coordinates of all the objects

```R
data(spARI_example_data)
```

Run spARI function to compute spRI and spARI. The meaning of each argument in the main function is listed below.

* r_labels: the reference partition
* c_labels: the clustering partition
* coords: spatial coordinates of all the objects
* alpha_val: coefficient belongs to the open interval $(0, 1)$ to keep a positive gap between the maximal weight of the disagreement pair and the weight one of the agreement pair. Default is 0.8.
* print_time: Logical; if TRUE, the total execution time is printed. Default is FALSE.

```R
res_value1 = spARI(true_labels, c1_labels, coords)
res_value2 = spARI(true_labels, c2_labels, coords)
```

Users can simply run the code `example("spARI")` to carry out this example.

```R
library(BACT)
example("spARI")
# 1st method: spRI=0.926, spARI=0.672
# 2nd method: spRI=0.917, spARI=0.632
```



## Remarks

If you have any questions regarding this package, please contact Yinqiao Yan at [yinqiaoyan@bjut.edu.cn](mailto:yinqiaoyan@bjut.edu.cn).

