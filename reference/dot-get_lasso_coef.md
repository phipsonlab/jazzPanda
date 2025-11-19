# help function to get lasso coefficient for every cluster for a given model

help function to get lasso coefficient for every cluster for a given
model

## Usage

``` r
.get_lasso_coef(
  i_gene,
  gene_mt,
  vec_cluster,
  cluster_names,
  n_fold = 10,
  n_samples,
  sample_names
)
```

## Arguments

- i_gene:

  Name of the current gene

- gene_mt:

  A matrix contains the transcript count in each grid. Each row refers
  to a grid, and each column refers to a gene. The column names must be
  specified and refer to the genes. This can be the output from the
  function
  [`get_vectors`](https://phipsonlab.github.io/jazzPanda/reference/get_vectors.md).

- vec_cluster:

  A matrix of the spatial vectors for clusters.

- cluster_names:

  A vector of strings giving the name of clusters

- n_fold:

  Optional. A positive number giving the number of folds used for cross
  validation. This parameter will pass to
  [`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html) to
  calculate a penalty term for every gene.

- n_samples:

  A positive number giving the number samples

- sample_names:

  A vector specifying the names for the sample

## Value

a list of two matrices with the following components

- `coef_df` :

  A matrix giving the lasso coefficient of each cluster

- `lambda.1se` :

  the lambda.1se value of best fitted model
