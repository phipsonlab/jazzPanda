# A small section of Xenium human breast cancer rep2.

A SpatialExperiment object containing the coordinates for every negative
control detection for rep2_sub

## Usage

``` r
data(rep2_sub)
```

## Format

A SpatialExperiment object with 20 genes and 1829 cells. The molecules
assay slot is a BumpyDataFrameMatrix obejct. Can retrieve DataFrame
version by calling
\`BumpyMatrix::unsplitAsDataFrame(molecules(rep2_sub))\`. The molecules
slot contains:

- x:

  x coordinates

- y:

  y coordiantes

- feature_name:

  transcript name

## Source

\<https://cf.10xgenomics.com/samples/xenium/1.0.1/
Xenium_FFPE_Human_Breast_Cancer_Rep2/
Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip\>

## Value

SpatialExperiment
