% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_data.R
\name{convert_data}
\alias{convert_data}
\title{Convert SingleCellExperiment/SpatialExperiment/SpatialFeatureExperiment 
objects to list object for jazzPanda.}
\usage{
convert_data(x, sample_names)
}
\arguments{
\item{x}{a SingleCellExperiment or SpatialExperiment or 
SpatialFeatureExperiment object}

\item{sample_names}{a vector of strings giving the sample names}
}
\value{
outputs a list object with the following components
\item{trans_lst }{ A list of named dataframes. Each dataframe refers to 
one sample and shows the transcript detection coordinates for each gene. 
The name matches the input sample_names}
\item{cm_lst }{A list of named dataframes containing the count matrix for
each sample. The name matches the input sample_names}
}
\description{
This function takes an object of class SingleCellExperiment, 
SpatialExperiment or SpatialFeatureExperimentreturns and returns a list
object that is expected for the \code{get_vector} functions.
}
\examples{

library(SingleCellExperiment)
library(SpatialExperiment)
# A SingleCellExperiment object
set.seed(200)
counts_sp1 <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
counts_sp2 <- matrix(rpois(100, lambda = 5), ncol=10, nrow=10)
sce <- SingleCellExperiment(list(sp1=counts_sp1, sp2=counts_sp2))
sce_example_lst <- convert_data(sce, sample_names = c("sp1","sp2"))

# A SpatialExperiment object
n <- 10
y <- matrix(rpois(200, lambda = 2),nrow = n, ncol = 2*n)
cd <- DataFrame(x = seq(2*n), y = seq(2*n))
spe1 <- SpatialExperiment(
    assays = list(counts = y),
    colData = cd, 
    sample_id="sample1",
    spatialCoordsNames = c("x", "y"))
se_example_lst <- convert_data(spe1, sample_names = c("sample1"))
## Multiple sample scenario
spe2 <- SpatialExperiment(
    assays = list(counts = matrix(rpois(200, lambda = 2),
                                  nrow = n, ncol = 2*n)),
    colData = cd, 
    sample_id="sample2",
    spatialCoordsNames = c("x", "y"))
se_example_lst <- convert_data(cbind(spe1, spe2), 
                               sample_names = c("sample1", "sample2"))

}
