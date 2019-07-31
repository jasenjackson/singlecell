
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
library(MASS)
library(Seurat)
```

Individual 4
------------

Import the count matrix.

``` r
Ind4.mtx <- read.table(file = "/projects/b1101/Jasen/data/Breast_scRNA_Kessenbrock_Ind4.txt", header = TRUE, row.names=1)
head(Ind4.mtx[,1:3])
```

    ##               Ind4_AAACATACGTACAC Ind4_AAACATTGCCTCCA Ind4_AAACATTGTGAAGA
    ## RP11-34P13.3                    0                   0                   0
    ## FAM138A                         0                   0                   0
    ## OR4F5                           0                   0                   0
    ## RP11-34P13.7                    0                   0                   0
    ## RP11-34P13.8                    0                   0                   0
    ## RP11-34P13.14                   0                   0                   0

Calculate statistics for each gene and cell. Plot QC graphs

``` r
Ind4.cellstats <- QCstats(Ind4.mtx, verbose=FALSE)
QCplot(Ind4.cellstats, hline=800, vline=18000) 
```

![](figures/QC_Kessenbrock-unnamed-chunk-4-1.png)

Load count matrix and cellstats metadata into a Seurat object. Perform QC filtering of cells and genes

``` r
Ind4.cellstats <- QCseurat(Ind4.cellstats) # formats cellstats matrix for Seurat
Ind4 <- CreateSeuratObject(count=Ind4.mtx, project="Breast_Ind4", meta.data=Ind4.cellstats)
Ind4 <- subset(x = Ind4, subset = percent.cg > 90 & gene_count > 800 & umi_count < 18000)
ncol(Ind4.mtx) - ncol(Ind4); ((ncol(Ind4.mtx)-ncol(Ind4)) / ncol(Ind4.mtx)) * 100
```

    ## [1] 151

    ## [1] 3.66861

151 cells (~3.7%) were filtered out. <br />

``` r
# save the seurat object and remove raw data from memory
saveRDS(Ind4, file="/projects/b1101/Jasen/data/BreastAtlas_Ind4.rds")
rm(Ind4.cellstats); rm(Ind4.mtx)
```
