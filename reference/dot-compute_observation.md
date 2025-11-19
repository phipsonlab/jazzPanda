# Compute observation statistic for permutation framework

Compute observation statistic for permutation framework

## Usage

``` r
.compute_observation(
  x,
  cluster_info,
  correlation_method,
  n_cores,
  test_genes,
  bin_type,
  bin_param,
  use_cm
)
```

## Arguments

- x:

  a named list (of transcript detection coordinates) or
  SingleCellExperiment or SpatialExperiment or SpatialFeatureExperiment
  object. If a named list is provided, every list element is a dataframe
  containing the transcript detection coordinates and column names must
  include "feature_name" (gene name), "x" (x coordinate), "y" (y
  coordinate). The list names must match samples in cluster_info.

- cluster_info:

  A dataframe/matrix containing the centroid coordinates, cluster label
  and sample for each cell.The column names must include "x" (x
  coordinate), "y" (y coordinate), "cluster" (cluster label) and
  "sample" (sample). It is strongly recommended to use syntactically
  valid names for columns clusters and samples. If invalid names are
  detected, the function
  [`make.names`](https://rdrr.io/r/base/make.names.html) will be
  employed to generate valid names. A message will also be displayed to
  indicate this change.

- correlation_method:

  A parameter pass to [`cor`](https://rdrr.io/r/stats/cor.html),
  indicating which correlation coefficient is to be computed. One of
  "pearson" (default), "kendall", or "spearman": can be abbreviated.

- n_cores:

  A positive number specifying number of cores used for parallelizing
  permutation testing. Default is one core (sequential processing).

- test_genes:

  A vector of strings giving the name of the genes you want to test
  correlation for.

- bin_type:

  A string indicating which bin shape is to be used for vectorization.
  One of "square" (default), "rectangle", or "hexagon".

- bin_param:

  A numeric vector indicating the size of the bin. If the `bin_type` is
  "square" or "rectangle", this will be a vector of length two giving
  the numbers of rectangular quadrats in the x and y directions. If the
  `bin_type` is "hexagonal", this will be a number giving the side
  length of hexagons. Positive numbers only.

- use_cm:

  A boolean value that specifies whether to create spatial vectors for
  genes using the count matrix and cell coordinates instead of the
  transcript coordinates when both types of information are available.
  The default setting is FALSE.

## Value

A named list with the following components

- `obs.stat` :

  A matrix contains the observation statistic for every gene and every
  cluster. Each row refers to a gene, and each column refers to a
  cluster

- `gene_mt` :

  contains the transcript count in each grid. Each row refers to a grid,
  and each column refers to a gene.
