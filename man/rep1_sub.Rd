% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jazzpanda-package.R
\docType{data}
\name{rep1_sub}
\alias{rep1_sub}
\title{A small section of Xenium human breast cancer rep1.}
\format{
A SpatialExperiment object with 20 genes and 1713 cells. 
The molecules assay slot is a BumpyDataFrameMatrix obejct. 
Can retrieve DataFrame version by calling
`BumpyMatrix::unsplitAsDataFrame(molecules(rep2_sub))`.
The molecules assay contains 79576 rows and 3 variables:
\describe{
    \item{x}{x coordinates}
    \item{y}{y coordiantes}
    \item{feature_name}{transcript name}
}
}
\source{
<https://cf.10xgenomics.com/samples/xenium/1.0.1/
Xenium_FFPE_Human_Breast_Cancer_Rep1/
Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip>
}
\usage{
data(rep1_sub)
}
\value{
SpatialExperiment
}
\description{
A SpatialExperiment object containing the coordinates for every transcript
}
\keyword{datasets}
