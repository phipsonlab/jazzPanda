#' A small section of Xenium human breast cancer rep1.
#'
#' A SpatialExperiment object containing the coordinates for every transcript
#'
#' @format A SpatialExperiment object with 20 genes and 1713 cells. 
#' The molecules assay slot is a BumpyDataFrameMatrix obejct. 
#' Can retrieve DataFrame version by calling
#' `BumpyMatrix::unsplitAsDataFrame(molecules(rep2_sub))`.
#' The molecules assay contains 79576 rows and 3 variables:
#' \describe{
#'     \item{x}{x coordinates}
#'     \item{y}{y coordiantes}
#'     \item{feature_name}{transcript name}
#' }
#' @usage data(rep1_sub)
#' @return SpatialExperiment
#' @source <https://cf.10xgenomics.com/samples/xenium/1.0.1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip>
"rep1_sub"


#' Rep1 selected cells
#'
#' A data frame file containing the coordinates and cluster label for 
#' each cell of the selected subset of rep1.
#'
#' @format A data frame with 1705 rows and 6 variables:
#' \describe{
#'     \item{anno}{the provided cell type annotation}
#'     \item{cluster}{cluster label}
#'     \item{x}{x coordinates}
#'     \item{y}{y coordiantes}
#'     \item{cells}{cell id}
#'     \item{sample}{sample id}
#' }
#' @usage data(rep1_clusters)
#' @return  data frame
#' @source <https://cf.10xgenomics.com/samples/xenium/1.0.1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip>
"rep1_clusters"

#' Rep1 negative control genes within the selected region.
#'
#' A SpatialExperiment object containing the coordinates for every negative 
#' control detection for rep1_sub  
#'
#' @format A SpatialExperiment object. 
#' The molecules assay slot is a BumpyDataFrameMatrix obejct. 
#' Can retrieve DataFrame version by calling
#' `BumpyMatrix::unsplitAsDataFrame(molecules(rep1_neg))`. 
#' The molecules slot contains: 
#' \describe{
#'     \item{x}{x coordinates}
#'     \item{y}{y coordiantes}
#'     \item{feature_name}{negative control probe name}
#'     \item{category}{negative control category}
#' }
#' @usage data(rep1_neg)
#' @return SpatialExperiment 
#' @source <https://cf.10xgenomics.com/samples/xenium/1.0.1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_
#' Cancer_Rep1_outs.zip>
"rep1_neg"


#' A small section of Xenium human breast cancer rep2.
#' 
#' A SpatialExperiment object containing the coordinates for every negative 
#' control detection for rep2_sub  
#' @format A SpatialExperiment object with 20 genes and 1829 cells.
#' The molecules assay slot is a BumpyDataFrameMatrix obejct. 
#' Can retrieve DataFrame version by calling
#' `BumpyMatrix::unsplitAsDataFrame(molecules(rep2_sub))`.
#' The molecules slot contains:
#' \describe{
#'     \item{x}{x coordinates}
#'     \item{y}{y coordiantes}
#'     \item{feature_name}{transcript name}
#' }
#' @usage data(rep2_sub)
#' @return SpatialExperiment
#' @source <https://cf.10xgenomics.com/samples/xenium/1.0.1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep2/
#' Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip>
"rep2_sub"


#' Rep2 selected cells
#'
#' A csv file containing the coordinates and cluster label for each cell of the
#' selected subset of rep2.
#'
#' @format A data frame with 1815 rows and and 6 variables:
#' \describe{
#'     \item{anno}{the provided cell type annotation}
#'     \item{cluster}{cluster label}
#'     \item{x}{x coordinates}
#'     \item{y}{y coordiantes}
#'     \item{cells}{cell id}
#'     \item{sample}{sample id}
#' }
#' @usage data(rep2_clusters)
#' @return  data frame
#' @source <https://cf.10xgenomics.com/samples/xenium/1.0.1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep2/
#' Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip>
"rep2_clusters"

#' Rep2 negative control genes within the selected region.
#'
#' A data frame containing the coordinates for every negative control detection 
#' for rep2
#' @format A SpatialExperiment object. 
#' The molecules assay slot is a BumpyDataFrameMatrix obejct. 
#' Can retrieve DataFrame version by calling
#' `BumpyMatrix::unsplitAsDataFrame(molecules(rep2_neg))`.
#' The molecules slot contains:
#' \describe{
#'     \item{x}{x coordinates}
#'     \item{y}{y coordiantes}
#'     \item{feature_name}{negative control probe name}
#'.    \item{category}{negative control category}
#' }
#' @usage data(rep2_neg)
#' @return SpatialExperiment
#' @source <https://cf.10xgenomics.com/samples/xenium/1.0.1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep2/
#' Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip>
"rep2_neg"
