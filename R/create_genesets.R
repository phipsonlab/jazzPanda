#' Convert the coordinates of set of genes into vectors.
#'
#' @param x a named list (of transcript detection coordinates) or 
#' SingleCellExperiment or SpatialExperiment or 
#' SpatialFeatureExperiment object. If a named list is provided, every list
#' element is a dataframe containing the transcript detection
#' coordinates and column names must include 
#' "feature_name" (nagative control name), "x" (x coordinate) and
#'  "y" (y coordinate). The list names must match samples in cluster_info. 
#' @param name_lst A named list of strings giving the name of features that are
#' treated as background.
#' @param cluster_info A dataframe/matrix containing the centroid coordinates,
#' cluster label and sample for each cell.The column names must include
#' "x" (x coordinate), "y" (y coordinate),
#' "cluster" (cluster label) and "sample" (sample). It is strongly recommended 
#' to use syntactically valid names for columns clusters and samples. 
#' If invalid names are detected, the function \code{\link{make.names}} will be 
#' employed to generate valid names. A message will also be displayed to 
#' indicate this change.
#' This parameter is only used when use_cm parameter is TURE
#' @param sample_names a vector of strings giving the sample names 
#' @param bin_type A string indicating which bin shape is to be used for
#' vectorization. One of "square" (default), "rectangle", or "hexagon".
#' @param bin_param A numeric vector indicating the size of the bin. If the
#' \code{bin_type} is "square" or "rectangle", this will be a vector of length
#' two giving the numbers of rectangular quadrats in the x and y directions. If
#' the \code{bin_type} is "hexagonal", this will be a number giving the side
#' length of hexagons. Positive numbers only.
#' @param use_cm A boolean value that specifies whether to create spatial 
#' vectors for genes using the count matrix and cell coordinates instead of 
#' the transcript coordinates when both types of information are available. 
#' The default setting is FALSE.
#' @param n_cores A positive number specifying number of cores used for
#' parallelizing permutation testing. Default is one core
#' (sequential processing).
#' @return a list of two matrices with the following components
#' \item{\code{gene_mt}  }{contains the transcript count in each grid.
#' Each row refers to a grid, and each column refers to a gene.}
#' \item{\code{cluster_mt}  }{contains the number of cells in a specific
#' cluster in each grid. Each row refers to a grid, and each column refers
#' to a cluster.}
#' The row order of \code{gene_mt} matches the row order of \code{cluster_mt}.

#' @param bin_type A string indicating which bin shape is to be used for
#' vectorization. One of "square" (default), "rectangle", or "hexagon".
#' @param bin_param A numeric vector indicating the size of the bin. If the
#' \code{bin_type} is "square" or "rectangle", this will be a vector of length
#' two giving the numbers of rectangular quadrats in the x and y directions. If
#' the \code{bin_type} is "hexagonal", this will be a number giving the side
#' length of hexagons. Positive numbers only.
#' @param cluster_info A dataframe/matrix containing the centroid coordinates,
#' cluster and sample label for each cell.The column names must include
#' "x" (x coordinate), "y" (y coordinate),
#' "cluster" (cluster label) and "sample" (sample).
#'
#' @return A matrix contains the sum count in each grid.
#' Each row refers to a grid, each column refers to a set in \code{name_lst}.
#' The column name will match the names in \code{name_lst}.
#' @export
#'
#' @examples
#' library(SpatialExperiment)
#' set.seed(15)
#' trans = as.data.frame(rbind(cbind(x = runif(10, min=1, max=10),
#'                          y = runif(10, min=1, max=10),
#'                                 feature_name="A"),
#'                          cbind(x = runif(5, min=10, max=24),
#'                                y = runif(5, min=1, max=10),
#'                                feature_name="B"),
#'                          cbind(x = runif(10, min=10, max=24),
#'                                y = runif(10, min=10, max=24),
#'                                feature_name="C")))
#' trans$x = as.numeric(trans$x)
#' trans$y = as.numeric(trans$y)
#' trans$cell = sample(c("cell1","cell2","cell2"),replace=TRUE,
#'                         size=nrow(trans))
#' # create SpatialExperiment object
#' trans_mol <- BumpyMatrix::splitAsBumpyMatrix(
#'     trans[, c("x", "y")], 
#'     row = trans$feature_name, col = trans$cell )

#' rep1_spe<- SpatialExperiment(
#'      assays = list(molecules = trans_mol),sample_id ="sample1" )

#' geneset_res <- create_genesets(x=rep1_spe, sample=c("sample1"),
#'                              name_lst=list(dummy_A=c("A","C"),
#'                                              dummy_B=c("A","B","C")),
#'                              bin_type="square",
#'                              bin_param=c(2,2),cluster_info=NULL)
#'
create_genesets<-function(x, name_lst, 
                        cluster_info, sample_names, bin_type, bin_param,
                        use_cm = FALSE, n_cores=1){
    primary_class <- class(x)[1]
    if (primary_class == "list"){
        uni_genes <- unique(unlist(lapply(x, function(df) df$feature_name)))
        for (sp in sample_names){
            sub_x <- x[[sp]]
            req_cols <- c("feature_name","x","y")
            if (length(setdiff(req_cols, colnames(sub_x))) != 0){
                stop("Invalid columns in input x. Every list element in x must
            contain columns 'feature_name','x', 'y', 
                for every transcipt")} }
    if (length(setdiff(as.vector(unique(unlist(name_lst))),uni_genes))>0){
            stop("Invalid name_lst, 
                can not find genes in the input name_lst from x")}
    test_genes <- intersect(as.vector(unique(unlist(name_lst))), uni_genes)
    }else {
    if (length(setdiff(as.vector(unique(unlist(name_lst))), row.names(x)))>0){
            stop("Invalid name_lst, can not find genes in the 
                input name_lst from x")
    }
    test_genes <- intersect(as.vector(unique(unlist(name_lst))), row.names(x))
    }
    
    res_lst <- get_vectors(x=x, cluster_info=cluster_info, 
                sample_names=sample_names, bin_type=bin_type, 
                bin_param=bin_param,test_genes=test_genes, 
                use_cm = use_cm, n_cores=n_cores)
    gene_vecs <- as.matrix(res_lst$gene_mt)
    
    vec_background <- as.data.frame(matrix(0, ncol=length(name_lst),
                            nrow=nrow(gene_vecs)))
    colnames(vec_background) <- names(name_lst)

    # iterate over each set over the name_lst
    for (nm in names(name_lst)){
        gs_genes <- unique(name_lst[[nm]])
        vec_background[,nm] <- rowSums(gene_vecs[,gs_genes, drop=FALSE])
    }
    return (vec_background)
}
