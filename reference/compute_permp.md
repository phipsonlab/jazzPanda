# Calculate a p-value for correlation with permutation.

This function will run permutation framework to compute a p-value for
the correlation between the vectorised genes and clusters each cluster
for one sample.

## Usage

``` r
compute_permp(
  x,
  cluster_info,
  perm.size,
  bin_type,
  bin_param,
  test_genes,
  correlation_method = "pearson",
  n_cores = 1,
  correction_method = "BH",
  use_cm = FALSE
)
```

## Arguments

- x:

  a named list (of transcript detection coordinates) or named
  SingleCellExperiment or named SpatialExperiment or named
  SpatialFeatureExperiment object. If a named list is provided, the list
  element is a dataframe containing the transcript detection coordinates
  and column names must include "feature_name" (gene name), "x" (x
  coordinate), "y" (y coordinate). The list name must match samples in
  cluster_info.

- cluster_info:

  A dataframe/matrix containing the centroid coordinates and cluster
  label for each cell.The column names should include "x" (x
  coordinate), "y" (y coordinate), and "cluster" (cluster label).

- perm.size:

  A positive number specifying permutation times

- bin_type:

  A string indicating which bin shape is to be used for vectorization.
  One of "square" (default), "rectangle", or "hexagon".

- bin_param:

  A numeric vector indicating the size of the bin. If the `bin_type` is
  "square" or "rectangle", this will be a vector of length two giving
  the numbers of rectangular quadrats in the x and y directions. If the
  `bin_type` is "hexagonal", this will be a number giving the side
  length of hexagons. Positive numbers only.

- test_genes:

  A vector of strings giving the name of the genes you want to test
  correlation for. `gene_mt`.

- correlation_method:

  A parameter pass to [`cor`](https://rdrr.io/r/stats/cor.html)
  indicating which correlation coefficient is to be computed. One of
  "pearson" (default), "kendall", or "spearman": can be abbreviated.

- n_cores:

  A positive number specifying number of cores used for parallelizing
  permutation testing. Default is one core (sequential processing).

- correction_method:

  A character string pass to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) specifying the
  correction method for multiple testing .

- use_cm:

  A boolean value that specifies whether to create spatial vectors for
  genes using the count matrix and cell coordinates instead of the
  transcript coordinates when both types of information are available.
  The default setting is FALSE.

## Value

An object of class 'cor_mg_result'. To access specific components of the
returned object:

- Use
  [`get_cor`](https://phipsonlab.github.io/jazzPanda/reference/get_cor.md)
  to retrieve the matrix of observed correlation coefficients.

- Use
  [`get_perm_p`](https://phipsonlab.github.io/jazzPanda/reference/get_perm_p.md)
  to access the matrix of raw permutation p-values.

- Use
  [`get_perm_adjp`](https://phipsonlab.github.io/jazzPanda/reference/get_perm_adjp.md)
  to obtain the matrix of adjusted permutation p-values.

## Details

To get a permutation p-value for the correlation between a gene and a
cluster, this function will permute the cluster label for each cell
randomly, and calculate correlation between the genes and permuted
clusters. This process will be repeated for `perm.size` times, and
permutation p-value is calculated as the probability of permuted
correlations larger than the observation correlation.

## Examples

``` r
library(SpatialExperiment)
#> Loading required package: SingleCellExperiment
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
library(BumpyMatrix)
set.seed(100)
# simulate coordinates for clusters
df_clA <- data.frame(x = rnorm(n=10, mean=20, sd=5),
                    y = rnorm(n=10, mean=20, sd=5), cluster="A")
df_clB <- data.frame(x = rnorm(n=10, mean=100, sd=5),
                    y = rnorm(n=10, mean=100, sd=5), cluster="B")
clusters <- rbind(df_clA, df_clB)
clusters$sample="sample1"
# simulate coordinates for genes
trans_info <- data.frame(rbind(cbind(x = rnorm(n=10, mean=20, sd=5),
                                    y = rnorm(n=10, mean=20, sd=5),
                                    feature_name="gene_A1"),
                    cbind(x = rnorm(n=10, mean=20, sd=5),
                                    y = rnorm(n=10, mean=20, sd=5),
                                    feature_name="gene_A2"),
                    cbind(x = rnorm(n=10, mean=100, sd=5),
                                    y = rnorm(n=10, mean=100, sd=5),
                                    feature_name="gene_B1"),
                    cbind(x = rnorm(n=10, mean=100, sd=5),
                                    y = rnorm(n=10, mean=100, sd=5),
                                    feature_name="gene_B2")))
trans_info$x<-as.numeric(trans_info$x)
trans_info$y<-as.numeric(trans_info$y)
trans_info$cell =  rep(paste("cell",1:20, sep=""), times=2)
mol <- BumpyMatrix::splitAsBumpyMatrix(
     trans_info[, c("x", "y")], 
     row = trans_info$feature_name, col = trans_info$cell )
spe_sample1 <- SpatialExperiment(
        assays = list(molecules = mol),sample_id ="sample1" )
set.seed(100)
corr_res <- compute_permp(x=spe_sample1,
             cluster_info=clusters,
             perm.size=10,
             bin_type="square",
             bin_param=c(2,2),
             test_genes=unique(trans_info$feature_name),
             correlation_method = "pearson",
             n_cores=1,
             correction_method="BH")
#> Correlation Method = pearson
#> Running 10 permutation in sequential
             
# raw permutation p-value
perm_p <- get_perm_p(corr_res)
# adjusted permutation p-value
adjusted_perm_p <- get_perm_adjp(corr_res)
# observed correlation 
obs_corr <- get_cor(corr_res)
```
