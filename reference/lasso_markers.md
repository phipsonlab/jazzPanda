# Find marker genes with spatial coordinates

This function will find the most spatially relevant cluster label for
each gene.

## Usage

``` r
lasso_markers(
  gene_mt,
  cluster_mt,
  sample_names,
  keep_positive = TRUE,
  background = NULL,
  n_fold = 10
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

- keep_positive:

  A logical flag indicating whether to return positively correlated
  clusters or not.

- background:

  Optional. A matrix providing the background information. Each row
  refers to a grid, and each column refers to one category of background
  information. Number of rows must equal to the number of rows in
  `gene_mt` and `cluster_mt`. Can be obtained by only providing
  coordinates matrices `cluster_info`. to function `get_vectors`.

- n_fold:

  Optional. A positive number giving the number of folds used for cross
  validation. This parameter will pass to
  [`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html) to
  calculate a penalty term for every gene.

## Value

An object of class 'glm_mg_result' To access specific components of the
returned object:

- Use
  [`get_top_mg`](https://phipsonlab.github.io/jazzPanda/reference/get_top_mg.md)
  to retrieve the top result data frame

- Use
  [`get_full_mg`](https://phipsonlab.github.io/jazzPanda/reference/get_full_mg.md)
  to retrieve full result data frame

## Details

This function will take the converted gene and cluster vectors from
function
[`get_vectors`](https://phipsonlab.github.io/jazzPanda/reference/get_vectors.md),
and return the most relevant cluster label for each gene. If there are
multiple samples in the dataset, this function will find shared markers
across different samples by including additional sample vectors in the
input `cluster_mt`.

This function treats all input cluster vectors as features, and create a
penalized linear model for one gene vector with lasso regularization.
Clusters with non-zero coefficient will be selected, and these clusters
will be used to formulate a generalised linear model for this gene
vector.

- If the input `keep_positive` is TRUE, the clusters with positive
  coefficient and significant p-value will be saved in the output matrix
  `lasso_full_result`. The cluster with a positive coefficient and the
  minimum p-value will be regarded as the most relevant cluster to this
  gene and be saved in the output matrix `lasso_result`.

- If the input `keep_positive` is FALSE, the clusters with negative
  coefficient and significant p-value will be saved in the output matrix
  `lasso_full_result`. The cluster with a negative coefficient and the
  minimum p-value will be regarded as the most relevant cluster to this
  gene and be saved in the output matrix `lasso_result`.

If there is no clusters with significant p-value, the a string "NoSig"
will be returned for this gene.

The parameter `background` can be used to capture unwanted noise pattern
in the dataset. For example, we can include negative control genes as a
background cluster in the model. If the most relevant cluster selected
by one gene matches the background "clusters", we will return "NoSig"
for this gene.

## See also

[`get_vectors`](https://phipsonlab.github.io/jazzPanda/reference/get_vectors.md)

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
top_result <- get_top_mg(lasso_res, coef_cutoff=0.05)
# the full result
full_result <- get_full_mg(lasso_res, coef_cutoff=0.05)
```
