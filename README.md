## spARI

Current version: 1.0 

The R package **spARI** is designed to compute two novel clustering evaluation metrics—the **spatially aware Rand index (spRI)** and its adjusted version (**spARI**)—for assessing spatial transcriptomics clustering. Unlike the traditional Rand index (RI) and adjusted Rand index (ARI), spRI and spARI incorporate spatial distance information into clustering evaluation. When comparing two partitions, spRI assigns a weight to each disagreement pair—two objects in the same cluster of one partition but in different clusters of the other—based on their spatial distance, allowing for a more refined distinction between disagreement pairs. The spRI value ranges between zero and one, while the spARI value is less than one with an expected value of zero. Higher spARI values indicate greater clustering accuracy and a more compact spatial structure. 

To implement this package, users only need to input the two partitions and the corresponding spatial coordinates, and the main function `spARI` efficiently computes spRI and spARI in a rapid and user-friendly manner. Additionally, the main function includes a built-in normalization step for spatial coordinates, eliminating the need for manual preprocessing. The **spARI** package is compatible with Windows, Linux, and macOS and can be easily installed across these platforms.



## Prerequisites and Installation

1. R version >= 4.1.3.
2. CRAN package: stats (>=4.1.3)
3. Install the package **spARI**.

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

The following code shows an example (the third simulation study in the manuscript) that runs the main function `spARI` in our package.

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
* alpha_val: coefficient belongs to the open interval (0, 1) to keep a positive gap between the maximal weight of the disagreement pair and the weight one of the agreement pair. Default is 0.8.

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


