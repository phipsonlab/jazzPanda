% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_permp.R
\name{compute_observation}
\alias{compute_observation}
\title{Compute observation statistic for permutation framework}
\usage{
compute_observation(
  data,
  cluster_info,
  correlation_method = "pearson",
  all_genes = all_genes,
  bin_type,
  bin_param,
  w_x,
  w_y
)
}
\arguments{
\item{data}{A list of matrices containing the coordinates of transcripts.}

\item{cluster_info}{A dataframe/matrix containing the centroid
coordinates and cluster label for each cell.The column names should
include "x" (x coordinate), "y" (y coordinate),
and "cluster" (cluster label).}

\item{correlation_method}{A parameter pass to \code{\link{cor}},
indicating which correlation coefficient is to be computed.
One of "pearson" (default), "kendall", or "spearman": can be abbreviated.}

\item{all_genes}{A vector of strings giving the name of the genes you want to
test correlation for.}

\item{bin_type}{A string indicating which bin shape is to be used for
vectorization. One of "square" (default), "rectangle", or "hexagon".}

\item{bin_param}{A numeric vector indicating the size of the bin. If the
\code{bin_type} is "square" or "rectangle", this will be a vector of length
two giving the numbers of rectangular quadrats in the x and y directions. If
the \code{bin_type} is "hexagonal", this will be a number giving the side
length of hexagons. Positive numbers only.}

\item{w_x}{a numeric vector of length two specifying the x coordinate
limits of enclosing box.}

\item{w_y}{a numeric vector of length two specifying the y coordinate
limits of enclosing box.}
}
\value{
A named list with the following components
\item{\code{obs.stat}  }{ A matrix contains the observation statistic for
every gene and every cluster. Each row refers to a gene, and each column
refers to a cluster}
\item{\code{gene_mt}  }{ contains the transcript count in each grid.
Each row refers to a grid, and each column refers to a gene.}
}
\description{
Compute observation statistic for permutation framework
}
