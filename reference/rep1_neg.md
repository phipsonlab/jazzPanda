# Rep1 negative control genes within the selected region.

A SpatialExperiment object containing the coordinates for every negative
control detection for rep1_sub

## Usage

``` r
data(rep1_neg)
```

## Format

A SpatialExperiment object. The molecules assay slot is a
BumpyDataFrameMatrix obejct. Can retrieve DataFrame version by calling
\`BumpyMatrix::unsplitAsDataFrame(molecules(rep1_neg))\`. The molecules
slot contains:

- x:

  x coordinates

- y:

  y coordiantes

- feature_name:

  negative control probe name

- category:

  negative control category

## Source

\<https://cf.10xgenomics.com/samples/xenium/1.0.1/
Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast\_
Cancer_Rep1_outs.zip\>

## Value

SpatialExperiment
