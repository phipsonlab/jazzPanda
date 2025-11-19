# A small section of Xenium human breast cancer rep1.

A SpatialExperiment object containing the coordinates for every
transcript

## Usage

``` r
data(rep1_sub)
```

## Format

A SpatialExperiment object with 20 genes and 1713 cells. The molecules
assay slot is a BumpyDataFrameMatrix obejct. Can retrieve DataFrame
version by calling
\`BumpyMatrix::unsplitAsDataFrame(molecules(rep2_sub))\`. The molecules
assay contains 79576 rows and 3 variables:

- x:

  x coordinates

- y:

  y coordiantes

- feature_name:

  transcript name

## Source

\<https://cf.10xgenomics.com/samples/xenium/1.0.1/
Xenium_FFPE_Human_Breast_Cancer_Rep1/
Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip\>

## Value

SpatialExperiment
