% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jazzpanda-package.R
\docType{data}
\name{rep2_sub}
\alias{rep2_sub}
\title{A small section of Xenium human breast cancer rep2.}
\format{
A SpatialExperiment object with 20 genes and 1829 cells.
The molecules assay slot is a BumpyDataFrameMatrix obejct. 
Can retrieve DataFrame version by calling
`BumpyMatrix::unsplitAsDataFrame(molecules(rep2_sub))`.
The molecules slot contains:
\describe{
    \item{x}{x coordinates}
    \item{y}{y coordiantes}
    \item{feature_name}{transcript name}
}
}
\source{
<https://cf.10xgenomics.com/samples/xenium/1.0.1/
Xenium_FFPE_Human_Breast_Cancer_Rep2/
Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip>
}
\usage{
data(rep2_sub)
}
\value{
SpatialExperiment
}
\description{
A SpatialExperiment object containing the coordinates for every negative 
control detection for rep2_sub
}
\keyword{datasets}
