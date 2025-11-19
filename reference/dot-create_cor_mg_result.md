# Create a marker gene result object for correlation approach

This function creates a structured output object named 'cor_mg_result'
for storing the permutation results. The object contains three matrices:

## Usage

``` r
.create_cor_mg_result(obs.stat, perm.pval, perm.pval.adj)
```

## Arguments

- obs.stat:

  A matrix containing the correlation coefficients for each pair of
  genes and cluster vectors.

- perm.pval:

  A matrix containing the raw permutation p-value for each pair of genes
  and cluster.

- perm.pval.adj:

  A matrix containing the adjusted permutation p-value for each pair of
  genes and cluster.

## Value

An S3 object of class 'cor_mg_result' which includes three matrices.
