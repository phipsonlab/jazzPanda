# Get top lasso result from glm_mg_result

Accessor function to retrieve the 'top_result' dataframe from an
'glm_mg_result' object.

## Usage

``` r
get_top_mg(obj, coef_cutoff = 0.05)
```

## Arguments

- obj:

  An 'glm_mg_result' object.

- coef_cutoff:

  A positive number giving the coefficient cutoff value. Genes whose top
  cluster showing a coefficient value smaller than the cutoff will be
  marked as non-marker genes ("NoSig"). Default is 0.05.

## Value

A data frame with detailed information for each gene and the most
relevant cluster label.

- `gene` Gene name

- `top_cluster` The name of the most relevant cluster after thresholding
  the coefficients.

- `glm_coef` The coefficient of the selected cluster in the generalised
  linear model.

- `pearson` Pearson correlation between the gene vector and the selected
  cluster vector.

- `max_gg_corr` A number showing the maximum pearson correlation for
  this gene vector and all other gene vectors in the input `gene_mt`

- `max_gc_corr` A number showing the maximum pearson correlation for
  this gene vector and every cluster vectors in the input `cluster_mt`

## Examples

``` r
library(SpatialExperiment)
set.seed(100)
#  simulate coordinates for clusters
df_clA <- data.frame(x = rnorm(n=100, mean=20, sd=5),
                 y = rnorm(n=100, mean=20, sd=5), cluster="A")
df_clB <- data.frame(x = rnorm(n=100, mean=100, sd=5),
                y = rnorm(n=100, mean=100, sd=5), cluster="B")

clusters <- rbind(df_clA, df_clB)
clusters$sample<-"sample1"

# simulate coordinates for genes
trans_info <- data.frame(rbind(cbind(x = rnorm(n=100, mean=20,sd=5),
                                y = rnorm(n=100, mean=20, sd=5),
                                 feature_name="gene_A1"),
                           cbind(x = rnorm(n=100, mean=20, sd=5),
                                 y = rnorm(n=100, mean=20, sd=5),
                                 feature_name="gene_A2"),
                           cbind(x = rnorm(n=100, mean=100, sd=5),
                                 y = rnorm(n=100, mean=100, sd=5),
                                 feature_name="gene_B1"),
                           cbind(x = rnorm(n=100, mean=100, sd=5),
                                 y = rnorm(n=100, mean=100, sd=5),
                                 feature_name="gene_B2")))
trans_info$x<-as.numeric(trans_info$x)
trans_info$y<-as.numeric(trans_info$y)
trans_info$cell<-sample(c("cell1","cell2","cell2"),replace=TRUE,
                        size=nrow(trans_info))
trans_mol <- BumpyMatrix::splitAsBumpyMatrix(
    trans_info[, c("x", "y")], 
    row = trans_info$feature_name, col = trans_info$cell )
spe<- SpatialExperiment(
     assays = list(molecules = trans_mol),sample_id ="sample1" )
vecs_lst <- get_vectors(x=spe,sample_names=c("sample1"),
                    cluster_info = clusters,
                    bin_type = "square",
                    bin_param = c(20,20),
                    test_genes =c("gene_A1","gene_A2","gene_B1","gene_B2"))
#> Warning: 1 point was rejected as lying outside the specified window
lasso_res <- lasso_markers(gene_mt=vecs_lst$gene_mt,
                        cluster_mt = vecs_lst$cluster_mt,
                        sample_names=c("sample1"),
                        keep_positive=TRUE,
                        background=NULL)
# the top result
top_result<- get_top_mg(lasso_res, coef_cutoff=0.05)
```
