# Convert the coordinates of set of genes into vectors.

Convert the coordinates of set of genes into vectors.

## Usage

``` r
create_genesets(
  x,
  name_lst,
  cluster_info,
  sample_names,
  bin_type,
  bin_param,
  use_cm = FALSE,
  n_cores = 1
)
```

## Arguments

- x:

  a named list (of transcript detection coordinates) or
  SingleCellExperiment or SpatialExperiment or SpatialFeatureExperiment
  object. If a named list is provided, every list element is a dataframe
  containing the transcript detection coordinates and column names must
  include "feature_name" (nagative control name), "x" (x coordinate) and
  "y" (y coordinate). The list names must match samples in cluster_info.

- name_lst:

  A named list of strings giving the name of features that are treated
  as background.

- cluster_info:

  A dataframe/matrix containing the centroid coordinates, cluster and
  sample label for each cell.The column names must include "x" (x
  coordinate), "y" (y coordinate), "cluster" (cluster label) and
  "sample" (sample).

- sample_names:

  a vector of strings giving the sample names

- bin_type:

  A string indicating which bin shape is to be used for vectorization.
  One of "square" (default), "rectangle", or "hexagon".

- bin_param:

  A numeric vector indicating the size of the bin. If the `bin_type` is
  "square" or "rectangle", this will be a vector of length two giving
  the numbers of rectangular quadrats in the x and y directions. If the
  `bin_type` is "hexagonal", this will be a number giving the side
  length of hexagons. Positive numbers only.

- use_cm:

  A boolean value that specifies whether to create spatial vectors for
  genes using the count matrix and cell coordinates instead of the
  transcript coordinates when both types of information are available.
  The default setting is FALSE.

- n_cores:

  A positive number specifying number of cores used for parallelizing
  permutation testing. Default is one core (sequential processing).

## Value

a list of two matrices with the following components

- `gene_mt` :

  contains the transcript count in each grid. Each row refers to a grid,
  and each column refers to a gene.

- `cluster_mt` :

  contains the number of cells in a specific cluster in each grid. Each
  row refers to a grid, and each column refers to a cluster.

The row order of `gene_mt` matches the row order of `cluster_mt`.

A matrix contains the sum count in each grid. Each row refers to a grid,
each column refers to a set in `name_lst`. The column name will match
the names in `name_lst`.

## Examples

``` r
library(SpatialExperiment)
set.seed(15)
trans = as.data.frame(rbind(cbind(x = runif(10, min=1, max=10),
                         y = runif(10, min=1, max=10),
                                feature_name="A"),
                         cbind(x = runif(5, min=10, max=24),
                               y = runif(5, min=1, max=10),
                               feature_name="B"),
                         cbind(x = runif(10, min=10, max=24),
                               y = runif(10, min=10, max=24),
                               feature_name="C")))
trans$x = as.numeric(trans$x)
trans$y = as.numeric(trans$y)
trans$cell = sample(c("cell1","cell2","cell2"),replace=TRUE,
                        size=nrow(trans))
# create SpatialExperiment object
trans_mol <- BumpyMatrix::splitAsBumpyMatrix(
    trans[, c("x", "y")], 
    row = trans$feature_name, col = trans$cell )
rep1_spe<- SpatialExperiment(
     assays = list(molecules = trans_mol),sample_id ="sample1" )
geneset_res <- create_genesets(x=rep1_spe, sample=c("sample1"),
                             name_lst=list(dummy_A=c("A","C"),
                                             dummy_B=c("A","B","C")),
                             bin_type="square",
                             bin_param=c(2,2),cluster_info=NULL)
```
