% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jazzpanda-package.R
\docType{data}
\name{rep2_neg}
\alias{rep2_neg}
\title{Rep2 negative control genes within the selected region.}
\format{
A SpatialExperiment object. 
The molecules assay slot is a BumpyDataFrameMatrix obejct. 
Can retrieve DataFrame version by calling
`BumpyMatrix::unsplitAsDataFrame(molecules(rep2_neg))`.
The molecules slot contains:
\describe{
    \item{x}{x coordinates}
    \item{y}{y coordiantes}
    \item{feature_name}{negative control probe name}
.    \item{category}{negative control category}
}
}
\source{
<https://cf.10xgenomics.com/samples/xenium/1.0.1/
Xenium_FFPE_Human_Breast_Cancer_Rep2/
Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip>
}
\usage{
data(rep2_neg)
}
\value{
SpatialExperiment
}
\description{
A data frame containing the coordinates for every negative control detection 
for rep2
}
\keyword{datasets}
