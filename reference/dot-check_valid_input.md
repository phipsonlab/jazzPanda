# helper function to check the inputs passed to marker detection function

helper function to check the inputs passed to marker detection function

## Usage

``` r
.check_valid_input(
  gene_mt,
  cluster_mt,
  sample_names,
  n_fold = 10,
  background = NULL
)
```

## Arguments

- gene_mt:

  A matrix contains the transcript count in each grid. Each row refers
  to a grid, and each column refers to a gene. The column names must be
  specified and refer to the genes. This can be the output from the
  function
  [`get_vectors`](https://phipsonlab.github.io/jazzPanda/reference/get_vectors.md).

- cluster_mt:

  A matrix contains the number of cells in a specific cluster in each
  grid. Each row refers to a grid, and each column refers to a cluster.
  The column names must be specified and refer to the clusters. Please
  do not assign integers as column names. This can be the output from
  the function
  [`get_vectors`](https://phipsonlab.github.io/jazzPanda/reference/get_vectors.md).

- sample_names:

  A vector specifying the names for the samples.

- n_fold:

  Optional. A positive number giving the number of folds used for cross
  validation. This parameter will pass to
  [`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html) to
  calculate a penalty term for every gene.

- background:

  Optional. A matrix providing the background information. Each row
  refers to a grid, and each column refers to one category of background
  information. Number of rows must equal to the number of rows in
  `gene_mt` and `cluster_mt`. Can be obtained by only providing
  coordinates matrices `cluster_info`. to function `get_vectors`.

## Value

a list of two matrices with the following components

- `n_clusters` :

  Number of clusters

- `cluster_names` :

  a vector of strings giving the name of the clusters
