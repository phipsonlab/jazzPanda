# Get permutation adjusted p value from cor_mg_result

Accessor function to retrieve the permutation adjusted p-value from an
'cor_mg_result' object.

## Usage

``` r
get_perm_adjp(obj)
```

## Arguments

- obj:

  An 'cor_mg_result' object.

## Value

A matrix contains the adjusted permutation p-value. Each row refers to a
gene, and each column refers to a cluster.

## Examples

``` r
library(SpatialExperiment)
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
# adjusted permutation p-value
adjusted_perm_p <- get_perm_adjp(corr_res)
```
