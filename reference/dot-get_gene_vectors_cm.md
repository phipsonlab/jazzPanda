# Create spatial vectors for genes from count matrix and cell coordinates

This function will build gene vectors with count matrix and cell
locations

## Usage

``` r
.get_gene_vectors_cm(
  cluster_info,
  cm_lst,
  bin_type,
  bin_param,
  test_genes,
  range_list
)
```

## Arguments

- cluster_info:

  A dataframe/matrix containing the centroid coordinates, cluster label
  and sample for each cell.The column names must include "x" (x
  coordinate), "y" (y coordinate), "cluster" (cluster label) and
  "sample" (sample).

- cm_lst:

  A list of named matrices containing the count matrix for each sample
  The name must match the sample column in `cluster_info`. If this input
  is provided, the `cluster_info` must be specified and contain an
  additional column "cell_id" to link cell location and count matrix.
  Default is NULL.

- bin_type:

  A string indicating which bin shape is to be used for vectorization.
  One of "square" (default), "rectangle", or "hexagon".

- bin_param:

  A numeric vector indicating the size of the bin. If the `bin_type` is
  "square" or "rectangle", this will be a vector of length two giving
  the numbers of rectangular quadrats in the x and y directions. If the
  `bin_type` is "hexagonal", this will be a number giving the side
  length of hexagons. Positive numbers only.

  For example:

  - `c(3, 4)` means 3 bins along the x-axis and 4 bins along the y-axis
    (a 3 × 4 grid).

  - `c(5, 5)` means 5 bins along the x-axis and 5 bins along the y-axis
    (a 5 × 5 grid).

- test_genes:

  A vector of strings giving the name of the genes you want to test.
  This will be used as column names for one of the result matrix
  `gene_mt`.

- range_list:

  A named list of spatial ranges for each sample. Each element should be
  a list with two components: `w_x` and `w_y`, which are numeric vectors
  of length 2 specifying the x- and y-axis ranges (e.g., from cell or
  transcript coordinates). The range is calculated with 5 within the
  window.

## Value

a matrix contains the transcript count in each grid. Each row refers to
a grid, and each column refers to a gene.
