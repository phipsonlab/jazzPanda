#' Convert SingleCellExperiment/SpatialExperiment/SpatialFeatureExperiment 
#' objects to list object for jazzPanda. 
#'
#' This function takes an object of class SingleCellExperiment, 
#' SpatialExperiment or SpatialFeatureExperimentreturns and returns a list
#' object that is expected for the \code{get_vector} functions. 
#' 
#' @param x a SingleCellExperiment or SpatialExperiment or 
#' SpatialFeatureExperiment object
#' @param sample_names a vector of strings giving the sample names 
#'
#' @return outputs a list object with the following components
#' \item{trans_lst }{ A list of named dataframes. Each dataframe refers to 
#' one sample and shows the transcript detection coordinates for each gene. 
#' The name matches the input sample_names}
#' \item{cm_lst }{A list of named dataframes containing the count matrix for
#' each sample. The name matches the input sample_names} 
#' 
#' @export
#' @importFrom methods is
#' @examples
#' 
#' library(SingleCellExperiment)
#' library(SpatialExperiment)
#' # A SingleCellExperiment object
#' set.seed(200)
#' counts_sp1 <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
#' counts_sp2 <- matrix(rpois(100, lambda = 5), ncol=10, nrow=10)
#' sce <- SingleCellExperiment(list(sp1=counts_sp1, sp2=counts_sp2))
#' sce_example_lst <- convert_data(sce, sample_names = c("sp1","sp2"))
#' 
#' # A SpatialExperiment object
#' n <- 10
#' y <- matrix(rpois(200, lambda = 2),nrow = n, ncol = 2*n)
#' cd <- DataFrame(x = seq(2*n), y = seq(2*n))
#' spe1 <- SpatialExperiment(
#'     assays = list(counts = y),
#'     colData = cd, 
#'     sample_id="sample1",
#'     spatialCoordsNames = c("x", "y"))
#' se_example_lst <- convert_data(spe1, sample_names = c("sample1"))
#' ## Multiple sample scenario
#' spe2 <- SpatialExperiment(
#'     assays = list(counts = matrix(rpois(200, lambda = 2),
#'                                   nrow = n, ncol = 2*n)),
#'     colData = cd, 
#'     sample_id="sample2",
#'     spatialCoordsNames = c("x", "y"))
#' se_example_lst <- convert_data(cbind(spe1, spe2), 
#'                                sample_names = c("sample1", "sample2"))
#' 
convert_data <- function(x, sample_names){
    n_samples <- length(sample_names)
    cm_lst <- vector("list", n_samples)
    names(cm_lst) <- sample_names
    trans_lst <- vector("list", n_samples)
    names(trans_lst) <- sample_names
    
    if (any(class(x) %in% c("matrix","array","list",
                            "character","data.frame"))){
        stop(sprintf("Unsupported input class. Please refer to the vignette 
            to simply create the required list"))
    }
    primary_class <- class(x)[1]
    if (primary_class == "SingleCellExperiment"){
    if (any(!sample_names %in% names(x@assays@data@listData))){
        stop("Input sample_names does not match the sample information in x.")}
        cm_lst <- x@assays@data@listData
        trans_lst <- NULL
    } else if (primary_class == "SpatialExperiment" | 
            primary_class ==  "SpatialFeatureExperiment"){
        if (!all(sample_names %in% unique(x$sample_id))){
            stop("Input sample_names does not match the 
                sample information in x.")}
        #all_cm <- x@assays@data$counts
        for (sp in sample_names){
            # x$imgData <- NULL
            sub_x <- x[,x$sample_id==sp]
            cm_lst[[sp]] <- sub_x@assays@data$counts }
            if (!is.null(x@assays@data$molecules)){
                for (sp in sample_names){
                    sub_x <- x[,x$sample_id==sp]
                    transcript_mlist <- sub_x@assays@data$molecules
                    transcript_nm<-row.names(sub_x@assays@data$molecules@proxy)
                    transcript_coords <- list()
                    for (nm in transcript_nm){
                    # make sure the transcript information contains columns: 
                    # feature_name, x, y
                        coords_df <- as.data.frame(unlist(transcript_mlist[nm]))
                        row.names(coords_df) <- seq_len(nrow(coords_df))
                        coords_df$feature_name <- nm
                        transcript_coords[[nm]] <- coords_df
                    }
                    tmp_data <- do.call("rbind", transcript_coords)
                    row.names(tmp_data) <- seq_len(nrow(tmp_data))
                    trans_lst[[sp]]<- tmp_data
                }
            }else{
                trans_lst <- NULL
            }
    }else{
        stop(sprintf("Unsupported input class. Please refer to the vignette 
                to simply create the required list")) }
    return (list(trans_lst = trans_lst, cm_lst = cm_lst))
}
