---
title: "QC Filtering of Breast Epithelial Atlas"
author: "Jasen Jackson"
date: "7/30/2019"
output:
  html_document:
    keep_md: TRUE
---




```r
# import packages and functions
setwd('/projects/b1101/Jasen')
source('singlecell/QCfunctions.R')

library(matrixStats)
library(scales)
library(ggplot2)
library(stringr)
library(dplyr)
```

## Import the count matrix
Import the count matrix for individual 4 from the Breast Epithelial atlas (Nguyen, et. al 2018; doi:[10.1038/s41467-018-04334-1](https://www.nature.com/articles/s41467-018-04334-1). The rows are genes and the columns are cells.


```r
print("hello world")
```

```
## [1] "hello world"
```
