# Create spatial vectors for clusters

Create spatial vectors for clusters

## Usage

``` r
.get_cluster_vectors(
  cluster_info,
  bin_length,
  bin_type,
  bin_param,
  range_list,
  sample_names
)
```

## Arguments

- cluster_info:

  A dataframe/matrix containing the centroid coordinates, cluster label
  and sample for each cell.The column names must include "x" (x
  coordinate), "y" (y coordinate), "cluster" (cluster label) and
  "sample" (sample).

- bin_length:

  A positive integer giving the length of total bins

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

- range_list:

  A named list of spatial ranges for each sample. Each element should be
  a list with two components: `w_x` and `w_y`, which are numeric vectors
  of length 2 specifying the x- and y-axis ranges (e.g., from cell or
  transcript coordinates). The range is calculated with 5 within the
  window.

- sample_names:

  a vector of strings giving the sample names

## Value

a matrix contains the cell count in each grid. Each row refers to a
grid, and each column refers to a cluster.
