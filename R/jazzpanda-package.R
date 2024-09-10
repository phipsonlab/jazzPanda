#' A small section of Xenium human breast cancer rep1.
#'
#' A list of named matrices
#'
#' @format A list of named matrices:
#' \describe{
#'     \item{cm}{the count matrix}
#'     \item{trans_info}{a data frame containing the coordinates for 
#'                      every transcript}
#' }
#' @usage data(rep1_sub)
#' @return List
#' @source \url{https://cf.10xgenomics.com/samples/xenium/1.0.1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip}
"rep1_sub"


#' Rep1 selected cells
#'
#' A csv file containing the coordinates and cluster label for each cell of the
#' selected subset of rep1.
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
#' @source \url{https://cf.10xgenomics.com/samples/xenium/1.0.1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip}
"rep1_clusters"

#' Rep1 negative control genes within the selected region.
#'
#' A list of named matrices
#'
#' @format A list of named matrices 
#' \describe{
#'     \item{cm}{the count matrix for negative control genes}
#'     \item{trans_info}{a data frame containing the coordinates 
#'                          for every transcript}
#'     \item{probe}{names of negative control probe genes}
#'     \item{codeword}{names of negative control codeword genes}
#' }
#' @usage data(rep1_neg)
#' @return List
#' @source \url{https://cf.10xgenomics.com/samples/xenium/1.0.1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_
#' Cancer_Rep1_outs.zip}
"rep1_neg"


#' A small section of Xenium human breast cancer rep2.
#'
#' A list of named matrices
#'
#' @format A list of named matrices 
#' \describe{
#'     \item{cm}{the count matrix}
#'     \item{trans_info}{a data frame containing the coordiantes 
#'                      for every transcript}
#' }
#' @usage data(rep2_sub)
#' @return List
#' @source \url{https://cf.10xgenomics.com/samples/xenium/1.0.1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep2/
#' Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip}
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
#' @source \url{https://cf.10xgenomics.com/samples/xenium/1.0.1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep2/
#' Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip}
"rep2_clusters"

#' Rep2 negative control genes within the selected region.
#'
#' A list of named matrices 
#'
#' @format A list of named matrices:
#' \describe{
#'     \item{cm}{the count matrix for negative control genes}
#'     \item{trans_info}{a data frame containing the coordinates for 
#'                      every transcript}
#'     \item{probe}{names of negative control probe genes}
#'     \item{codeword}{names of negative control codeword genes}
#' }
#' @usage data(rep2_neg)
#' @return List
#' @source \url{https://cf.10xgenomics.com/samples/xenium/1.0.1/
#' Xenium_FFPE_Human_Breast_Cancer_Rep2/
#' Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip}
"rep2_neg"
