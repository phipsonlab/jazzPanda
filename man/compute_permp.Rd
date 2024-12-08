% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_permp.R
\name{compute_permp}
\alias{compute_permp}
\title{Calculate a p-value for correlation with permutation.}
\usage{
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
  w_x,
  w_y,
  use_cm = FALSE
)
}
\arguments{
\item{x}{a SingleCellExperiment or SpatialExperiment or 
SpatialFeatureExperiment object}

\item{cluster_info}{A dataframe/matrix containing the centroid coordinates
and cluster label for each cell.The column names should include "x"
(x coordinate), "y" (y coordinate), and "cluster" (cluster label).}

\item{perm.size}{A positive number specifying permutation times}

\item{bin_type}{A string indicating which bin shape is to be used for
vectorization. One of "square" (default), "rectangle", or "hexagon".}

\item{bin_param}{A numeric vector indicating the size of the bin. If the
\code{bin_type} is "square" or "rectangle", this will be a vector of length
two giving the numbers of rectangular quadrats in the x and y directions. If
the \code{bin_type} is "hexagonal", this will be a number giving the side
length of hexagons. Positive numbers only.}

\item{test_genes}{A vector of strings giving the name of the genes you
want to test correlation for.
\code{gene_mt}.}

\item{correlation_method}{A parameter pass to \code{\link{cor}} indicating
which correlation coefficient is to be computed.
One of "pearson" (default), "kendall", or "spearman": can be abbreviated.}

\item{n_cores}{A positive number specifying number of cores used for
parallelizing permutation testing. Default is one core
(sequential processing).}

\item{correction_method}{A character string pass to \code{\link{p.adjust}}
specifying the correction method for multiple testing .}

\item{w_x}{a numeric vector of length two specifying the x coordinate
limits of enclosing box.}

\item{w_y}{a numeric vector of length two specifying the y coordinate
limits of enclosing box.}

\item{use_cm}{A boolean value that specifies whether to create spatial 
vectors for genes using the count matrix and cell coordinates instead of 
the transcript coordinates when both types of information are available. 
The default setting is FALSE.}
}
\value{
A named list with the following components
\item{\code{obs.stat}  }{ A matrix contains the observation statistic for
every gene and every cluster. Each row refers to a gene, and each column
refers to a cluster}
\item{\code{perm.arrays}  }{ A three dimensional array.
The first two dimensions represent the correlation between the genes and
permuted clusters. The third dimension refers to the different permutation
runs. }
\item{\code{perm.pval}  }{A matrix contains the raw permutation p-value.
Each row refers to a gene, and each column refers to a cluster}
\item{\code{perm.pval.adj}  }{A matrix contains the adjusted permutation
p-value. Each row refers to a gene, and each column refers to a cluster}
}
\description{
This function will run permutation framework to compute a p-value for the
correlation between the vectorised genes and clusters each cluster for one 
sample.
}
\details{
To get a permutation p-value for the correlation between a gene
and a cluster, this function will permute the cluster label for
each cell randomly, and calculate correlation between the genes and
permuted clusters. This process will be repeated for \code{perm.size}
times, and permutation p-value is calculated as the probability of
permuted correlations larger than the observation correlation.
}
\examples{
library(SpatialExperiment)
library(BumpyMatrix)
set.seed(100)
# simulate coordinates for clusters
df_clA = data.frame(x = rnorm(n=10, mean=20, sd=5),
                    y = rnorm(n=10, mean=20, sd=5), cluster="A")
df_clB = data.frame(x = rnorm(n=10, mean=100, sd=5),
                    y = rnorm(n=10, mean=100, sd=5), cluster="B")
clusters = rbind(df_clA, df_clB)
clusters$sample="sample1"
# simulate coordinates for genes
trans_info = data.frame(rbind(cbind(x = rnorm(n=10, mean=20, sd=5),
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
trans_info$x=as.numeric(trans_info$x)
trans_info$y=as.numeric(trans_info$y)
trans_info$cell =  rep(paste("cell",1:20, sep=""), times=2)
mol <- BumpyMatrix::splitAsBumpyMatrix(
     trans_info[, c("x", "y")], 
     row = trans_info$feature_name, col = trans_info$cell )
spe_sample1 <- SpatialExperiment(
        assays = list(molecules = mol),sample_id ="sample1" )
w_x =  c(min(floor(min(trans_info$x)),
             floor(min(clusters$x))),
         max(ceiling(max(trans_info$x)),
             ceiling(max(clusters$x))))
w_y =  c(min(floor(min(trans_info$y)),
             floor(min(clusters$y))),
         max(ceiling(max(trans_info$y)),
             ceiling(max(clusters$y))))
set.seed(100)
perm_p_lst = compute_permp(x=spe_sample1,
             cluster_info=clusters,
             perm.size=10,
             bin_type="square",
             bin_param=c(2,2),
             test_genes=unique(trans_info$feature_name),
             correlation_method = "pearson",
             n_cores=1,
             correction_method="BH",
             w_x=w_x ,
             w_y=w_y)
}
