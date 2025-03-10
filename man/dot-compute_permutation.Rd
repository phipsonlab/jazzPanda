% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_permp.R
\name{.compute_permutation}
\alias{.compute_permutation}
\title{Compute permutation statistics for permutation framework}
\usage{
.compute_permutation(
  cluster_info,
  perm.size = 1000,
  correlation_method = "pearson",
  bin_type,
  bin_param,
  n_cores = 1,
  w_x,
  w_y,
  gene_mt,
  cluster_names
)
}
\arguments{
\item{cluster_info}{A dataframe/matrix containing the centroid
coordinates and cluster label for each cell.The column names
should include "x" (x coordinate), "y" (y coordinate), and
"cluster" (cluster label).}

\item{perm.size}{A positive number specifying permutation times}

\item{correlation_method}{A parameter pass to \code{\link{cor}},
indicating which correlation coefficient is to be computed.
One of "pearson" (default), "kendall", or "spearman": can be abbreviated.}

\item{bin_type}{A string indicating which bin shape is to be used for
vectorization. One of "square" (default), "rectangle", or "hexagon".}

\item{bin_param}{A numeric vector indicating the size of the bin. If the
\code{bin_type} is "square" or "rectangle", this will be a vector of length
two giving the numbers of rectangular quadrats in the x and y directions. If
the \code{bin_type} is "hexagonal", this will be a number giving the side
length of hexagons. Positive numbers only.}

\item{n_cores}{A positive number specifying number of cores used for
parallelizing permutation testing. Default is one core
(sequential processing).}

\item{w_x}{a numeric vector of length two specifying the x coordinate
limits of enclosing box.}

\item{w_y}{a numeric vector of length two specifying the y coordinate
limits of enclosing box.}

\item{gene_mt}{A matrix contains the transcript count in each grid.
Each row refers to a grid, and each column refers to a gene.}

\item{cluster_names}{A list of strings giving the name and order of the
clusters}
}
\value{
A matrix with permutation statistics
}
\description{
Compute permutation statistics for permutation framework
}
