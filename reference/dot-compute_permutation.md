# Compute permutation statistics for permutation framework

Compute permutation statistics for permutation framework

## Usage

``` r
.compute_permutation(
  cluster_info,
  perm.size = 1000,
  correlation_method = "pearson",
  bin_type,
  bin_param,
  n_cores = 1,
  gene_mt,
  cluster_names,
  window_range
)
```

## Arguments

- cluster_info:

  A dataframe/matrix containing the centroid coordinates and cluster
  label for each cell.The column names should include "x" (x
  coordinate), "y" (y coordinate), and "cluster" (cluster label).

- perm.size:

  A positive number specifying permutation times

- correlation_method:

  A parameter pass to [`cor`](https://rdrr.io/r/stats/cor.html),
  indicating which correlation coefficient is to be computed. One of
  "pearson" (default), "kendall", or "spearman": can be abbreviated.

- bin_type:

  A string indicating which bin shape is to be used for vectorization.
  One of "square" (default), "rectangle", or "hexagon".

- bin_param:

  A numeric vector indicating the size of the bin. If the `bin_type` is
  "square" or "rectangle", this will be a vector of length two giving
  the numbers of rectangular quadrats in the x and y directions. If the
  `bin_type` is "hexagonal", this will be a number giving the side
  length of hexagons. Positive numbers only.

- n_cores:

  A positive number specifying number of cores used for parallelizing
  permutation testing. Default is one core (sequential processing).

- gene_mt:

  A matrix contains the transcript count in each grid. Each row refers
  to a grid, and each column refers to a gene.

- cluster_names:

  A list of strings giving the name and order of the clusters

- window_range:

  A list of spatial ranges for x and y. This list contains two
  components: `w_x` and `w_y`, which are numeric vectors of length 2
  specifying the x- and y-axis ranges (e.g., from cell or transcript
  coordinates).

## Value

A matrix with permutation statistics
