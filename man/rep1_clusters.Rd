% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jazzpanda-package.R
\docType{data}
\name{rep1_clusters}
\alias{rep1_clusters}
\title{Rep1 selected cells}
\format{
A data frame with 1705 rows and 6 variables:
\describe{
    \item{anno}{the provided cell type annotation}
    \item{cluster}{cluster label}
    \item{x}{x coordinates}
    \item{y}{y coordiantes}
    \item{cells}{cell id}
    \item{sample}{sample id}
}
}
\source{
<https://cf.10xgenomics.com/samples/xenium/1.0.1/
Xenium_FFPE_Human_Breast_Cancer_Rep1/
Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip>
}
\usage{
data(rep1_clusters)
}
\value{
data frame
}
\description{
A csv file containing the coordinates and cluster label for each cell of the
selected subset of rep1.
}
\keyword{datasets}
