
\name{jazzPanda-package}
\alias{jazzPanda-package}
\alias{jazzPanda}
\docType{package}
\title{
jazzPanda: A hybrid approach to find spatially relevant marker genes
in image-based spatial transcriptomics data
}
\description{
jazzPanda pacakge provides hybrid approaches to prioritize marker genes 
that uses the spatial coordinates of gene detections and cells making up 
clusters. 
We propose a binning approach \code{\link{get_vectors}}
that summarises the number of genes and cells
within a cluster as spatial vectors. We have developed two approaches to 
detect and prioritize marker genes. 
The first approach \code{\link{compute_permp}} is based on correlation 
coefficients between genes and cluster spatial vectors, where significance 
of the marker genes are assessed through permutation. 
The second approach \code{\link{lasso_markers}}
is based on lasso regularisation and linear modeling of our defined spatial 
vectors. This second approach is more flexible and can account for multiple 
samples and background noise. 
}

\author{
Melody Jin \email{jin.m@wehi.edu.au}
}

\keyword{package}
