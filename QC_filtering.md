
QC Filtering of Breast Epithelial Atlas
=======================================

``` r
# import packages and functions
setwd('/projects/b1101/Jasen')
source('singlecell/QCfunctions.R')

library(matrixStats)
library(scales)
library(ggplot2)
library(stringr)
library(dplyr)
```

Import the count matrix
-----------------------

Import the count matrix for iiindividual 4 from the Breast Epithelial atlas (Nguyen, et. al 2018; <doi:%5B10.1038/s41467-018-04334-1%5D(https://www.nature.com/articles/s41467-018-04334-1)>. <br /> <br /> The rows are genes and the columns are cells.

``` r
ggplot(midwest, aes(x=area, y=poptotal)) + geom_point()
```

![](images/qc-unnamed-chunk-3-1.png)
