[![DOI](https://zenodo.org/badge/724929670.svg)](https://zenodo.org/doi/10.5281/zenodo.10360070)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/phipsonlab/jazzPanda/blob/main/LICENSE)
[![R-CMD-check](https://github.com/phipsonlab/jazzPanda/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/phipsonlab/jazzPanda/actions/workflows/R-CMD-check.yaml)  

[![Logo](https://github.com/phipsonlab/jazzPanda/blob/main/inst/images/jazzPanda_logo.png)](https://github.com/phipsonlab/jazzPanda)

### 

jazzPanda: A hybrid approach to find spatially relevant marker genes

in image-based spatial transcriptomics data

  
[**Explore the vignette
»**](https://phipsonlab.github.io/jazzPanda/articles/jazzPanda.html)  
  

## Introduction

The jazzPanda package contains functions to find marker genes based on
the spatial coordinates for imaging-based spatial transcriptomics
technologies including 10x Xenium, NanoString CosMx and Vizgen MERSCOPE.

## Installation

To install jazzPanda from Bioconductor, use the following commands:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("jazzPanda")
```

To install jazzPanda from github, use the following commands:

``` r
library(devtools)
devtools::install_github("phipsonlab/jazzPanda")
```

In order to view the vignette for jazzPanda use the following command:

``` r
browseVignettes("jazzPanda")
```
