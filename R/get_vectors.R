#' Create spatial vectors for genes from transcript coordinates
#'
#' @description
#' This function will build gene vectors based on the transcript coordinates of
#' every gene
#' @param trans_lst If specified, it is a list of named dataframes. 
#' Each dataframe refers to one sample and shows the transcript detection 
#' coordinates for each gene. Optional parameter.
#' @param test_genes A vector of strings giving the name of the genes you 
#' want to test. This will be used as column names for one of the result matrix
#' \code{gene_mt}.
#' @param bin_type A string indicating which bin shape is to be used for
#' vectorization. One of "square" (default), "rectangle", or "hexagon".
#' @param bin_param A numeric vector indicating the size of the bin. If the
#' \code{bin_type} is "square" or "rectangle", this will be a vector of length
#' two giving the numbers of rectangular quadrats in the x and y directions. If
#' the \code{bin_type} is "hexagonal", this will be a number giving the side
#' length of hexagons. Positive numbers only.
#' @param bin_length A positive integer giving the length of total bins
#' @param w_x A numeric vector of length two specifying the x coordinate
#' limits of enclosing box.
#' @param w_y A numeric vector of length two specifying the y coordinate
#' limits of enclosing box.
#' @importFrom BiocParallel bplapply
#' @return a matrix contains the transcript count in each grid.
#' Each row refers to a grid, and each column refers to a gene.
get_gene_vectors_tr<- function(trans_lst, test_genes, bin_type, bin_param,
                                bin_length, w_x, w_y){
    n_genes <- length(test_genes)
    n_samples <- length(trans_lst)
    vec_gene_mt <- as.data.frame(matrix(0, ncol=n_genes,
                                        nrow=bin_length*n_samples))
    colnames(vec_gene_mt) <- test_genes
    if (bin_type == "hexagon"){
        w <- owin(xrange=w_x, yrange=w_y)
        H <- hextess(W=w, bin_param[1])
    }
    calculate_one_gene <- function(i_gene){
        vec_gene <- c()
        for (rpp in trans_lst){
            curr <- rpp[rpp$feature_name==i_gene,
                                    c("x","y")] %>% distinct()
            gene_ppp <- ppp(curr$x,curr$y,w_x, w_y)
            # create gene vector
            if (bin_type == "hexagon"){
                vec_g <- as.vector(t(quadratcount(gene_ppp, tess=H)))
            }else{
                vec_g <- as.vector(t(quadratcount(gene_ppp,
                                                bin_param[1],bin_param[2])))
            }
            vec_gene <- c(vec_gene, vec_g)
        }
        return(vec_gene)
    }

    # to calculate all genes in parallel
    result_lst <- bplapply(test_genes, calculate_one_gene)

    # convert the list of results to a data frame
    vec_gene_mt <- do.call(cbind, result_lst)
    
    # Convert to a data frame 
    vec_gene_mt <- as.data.frame(vec_gene_mt)
    colnames(vec_gene_mt) <- test_genes
    
    return (vec_gene_mt)
}


#' Create spatial vectors for genes from count matrix and cell coordinates
#'
#' @description
#' This function will build gene vectors with count matrix and cell locations
#' @param cluster_info A dataframe/matrix containing the centroid coordinates,
#' cluster label and sample for each cell.The column names must include
#' "x" (x coordinate), "y" (y coordinate),
#' "cluster" (cluster label) and "sample" (sample).
#' @param cm_lst A list of named matrices containing the count matrix for
#' each sample The name must match the sample column in \code{cluster_info}.
#' If this input is provided, the \code{cluster_info} must be specified and
#' contain an additional column "cell_id" to link cell location and
#' count matrix. Default is NULL.
#' @param bin_type A string indicating which bin shape is to be used for
#' vectorization. One of "square" (default), "rectangle", or "hexagon".
#' @param bin_param A numeric vector indicating the size of the bin. If the
#' \code{bin_type} is "square" or "rectangle", this will be a vector of length
#' two giving the numbers of rectangular quadrats in the x and y directions. If
#' the \code{bin_type} is "hexagonal", this will be a number giving the side
#' length of hexagons. Positive numbers only.
#' @param test_genes A vector of strings giving the name of the genes you 
#' want to test. This will be used as column names for one of the result matrix
#' \code{gene_mt}.
#' @param w_x A numeric vector of length two specifying the x coordinate
#' limits of enclosing box.
#' @param w_y A numeric vector of length two specifying the y coordinate
#' limits of enclosing box.
#' 
#' @importFrom BiocParallel bplapply
#'
#' @return a matrix contains the transcript count in each grid.
#' Each row refers to a grid, and each column refers to a gene.
#'
#' @importFrom spatstat.geom ppp
#' @importFrom spatstat.geom quadratcount
#' @importFrom spatstat.geom tileindex
#' @importFrom spatstat.geom owin
#' @importFrom spatstat.geom hextess
#' @importFrom dplyr distinct
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @importFrom spatstat.geom as.tess
#' @importFrom magrittr "%>%"
#'
get_gene_vectors_cm<- function(cluster_info, cm_lst, bin_type, bin_param,
                                test_genes, w_x, w_y){
    # binning
    bin_length <- 0
    if (bin_type == "hexagon"){
        w <- owin(xrange=w_x, yrange=w_y)
        H <- hextess(W=w, bin_param[1])
        bin_length <- length(H$tiles)
    }else{ bin_length <- bin_param[1] * bin_param[2] }
    n_genes <- length(test_genes)
    index_vec <- NULL
    # use count matrix to build gene vector matrix
    if (bin_type == "hexagon"){
        tile_indices <- tileindex(x=cluster_info$x, y=cluster_info$y, Z=H)
        cluster_info$index_vec  <- match(tile_indices,levels(tile_indices))
    }else{
        all_cell_intervals <- as.tess(quadratcount(ppp(cluster_info$x,
                        cluster_info$y, w_x, w_y),bin_param[1],bin_param[2]))
        cluster_info$index_x <- findInterval(cluster_info$x,
                                all_cell_intervals$xgrid, all.inside = TRUE, 
                                left.open=FALSE, rightmost.closed=TRUE )
        cluster_info$index_y <- findInterval(cluster_info$y,
                                all_cell_intervals$ygrid, all.inside = TRUE, 
                                left.open=FALSE, rightmost.closed=TRUE )
        cluster_info$index_y  <- bin_param[2] + 1 - cluster_info$index_y
        # column major
        ind_vec <-(cluster_info$index_y-1)* bin_param[2]+cluster_info$index_x
        cluster_info$index_vec <- ind_vec }
    n_samples <- length(cm_lst)
    # create gene vector
    vec_gene_mt<-as.data.frame(matrix(0,ncol=n_genes,nrow=bin_length*n_samples))
    colnames(vec_gene_mt) <- test_genes
    count_value <- NULL
    calculate_one_gene <- function(i_gene){
        vec_gene <- c()
        for (rp_nm in names(cm_lst)){
            cm <- cm_lst[[rp_nm]]
            vec_g <- rep(0, bin_length)
            i_gene_mt <- cluster_info[cluster_info$sample==rp_nm, ]
            i_gene_mt$count_value <- as.numeric(cm[i_gene,i_gene_mt$cell_id])
            i_gene_mt$index_vec <- factor(i_gene_mt$index_vec )
            added_count <- i_gene_mt %>% group_by(index_vec) %>%
                summarise(sum_ct = sum(count_value)) %>% data.frame
            idx_ve <- as.numeric(as.character((added_count$index_vec)))
            added_count$index_vec<-idx_ve
            vec_g[added_count[,"index_vec"]] <- added_count[,"sum_ct"]
            vec_gene <- c(vec_gene, vec_g)
        }
        return(vec_gene)
    }
    # calculate all genes in parallel
    result_lst <- bplapply(test_genes, calculate_one_gene)
    # Convert the list of results to a data frame
    vec_gene_mt <- do.call(cbind, result_lst)
    vec_gene_mt <- as.data.frame(vec_gene_mt)
    colnames(vec_gene_mt) <- test_genes
    return (vec_gene_mt)
}

#' Create spatial vectors for clusters
#'
#' @param cluster_info A dataframe/matrix containing the centroid coordinates,
#' cluster label and sample for each cell.The column names must include
#' "x" (x coordinate), "y" (y coordinate),
#' "cluster" (cluster label) and "sample" (sample).
#' @param bin_length A positive integer giving the length of total bins
#' @param bin_type A string indicating which bin shape is to be used for
#' vectorization. One of "square" (default), "rectangle", or "hexagon".
#' @param bin_param A numeric vector indicating the size of the bin. If the
#' \code{bin_type} is "square" or "rectangle", this will be a vector of length
#' two giving the numbers of rectangular quadrats in the x and y directions. If
#' the \code{bin_type} is "hexagonal", this will be a number giving the side
#' length of hexagons. Positive numbers only.
#' @param w_x A numeric vector of length two specifying the x coordinate
#' limits of enclosing box.
#' @param w_y A numeric vector of length two specifying the y coordinate
#' limits of enclosing box.
#' @param sample_names a vector of strings giving the sample names
#' @return a matrix contains the cell count in each grid.
#' Each row refers to a grid, and each column refers to a cluster.
#'
get_cluster_vectors<- function(cluster_info,bin_length,bin_type, bin_param,
                                w_x, w_y, sample_names){
    ava_sample_names <- unique(as.character(cluster_info$sample))
    if (length(setdiff(sample_names,ava_sample_names))>0){
        stop("Can not find input sample_names from input cluster_info")
    }
    n_clusters <- length(unique(cluster_info$cluster))
    n_samples <- length(sample_names)
    cluster_names <- unique(as.character(cluster_info$cluster))
    sample_cluster_nm <- paste(rep(sample_names, each=n_clusters),
                            cluster_names,sep="--")
    # create cluster dataframe
    vec_cluster <- as.data.frame(matrix(0, nrow=n_samples*bin_length,
                                        ncol=n_clusters))
    if (n_samples >1){
        sps <- as.data.frame(rep(sample_names, each=bin_length))
        colnames(sps) <- "samples"
        samples_variables <- model.matrix(~ samples- 1, data = sps)
        colnames(samples_variables) <- sample_names
        vec_cluster <- as.data.frame(cbind(vec_cluster,samples_variables))
        colnames(vec_cluster) <- c(as.character(cluster_names), sample_names)
    }else{
        vec_cluster <- as.data.frame(vec_cluster)
        colnames(vec_cluster) <- as.character(cluster_names)
    }
    # create cluster vector
    for (nm in sample_cluster_nm ){
        rp <- unlist(strsplit(nm, split="--"))[1]
        cl <- unlist(strsplit(nm, split="--"))[2]
        cluster_ppp <- ppp(cluster_info[cluster_info$sample==rp &
                                        cluster_info$cluster==cl,
                                        "x"],
                        cluster_info[cluster_info$sample==rp &
                                        cluster_info$cluster==cl,
                                        "y"], w_x, w_y)
        if (bin_type == "hexagon"){
            w <- owin(xrange=w_x, yrange=w_y)
            H <- hextess(W=w, bin_param[1])
            cm_cluster <- quadratcount(cluster_ppp, tess=H)
        }else{
            cm_cluster <- quadratcount(cluster_ppp, bin_param[1],bin_param[2])
        }
        if (n_samples >1){
            vec_cluster[vec_cluster[,rp]==1,cl] <- as.vector(t(cm_cluster))
        }else{
            vec_cluster[,cl] <- as.vector(t(cm_cluster))
        }
    }
    return (vec_cluster)
}


#' helper function to check the input of binning
#'
#' @param bin_type A string indicating which bin shape is to be used for
#' vectorization. One of "square" (default), "rectangle", or "hexagon".
#' @param bin_param A numeric vector indicating the size of the bin. If the
#' \code{bin_type} is "square" or "rectangle", this will be a vector of length
#' two giving the numbers of rectangular quadrats in the x and y directions. If
#' the \code{bin_type} is "hexagonal", this will be a number giving the side
#' length of hexagons. Positive numbers only.
#' @param w_x A numeric vector of length two specifying the x coordinate
#' limits of enclosing box.
#' @param w_y A numeric vector of length two specifying the y coordinate
#' limits of enclosing box.
#'
#' @return the length of total bins

check_binning<- function(bin_param, bin_type, w_x, w_y){
    # binning
    bin_length <- 0
    if (bin_type == "hexagon"){
        if (length(bin_param) != 1){
            stop("Invalid input bin_param, bin_param should be a vector
                of length 2 for hexagon bins")
        }
        w <- owin(xrange=w_x, yrange=w_y)
        H <- hextess(W=w, bin_param[1])
        bin_length <- length(H$tiles)
    }else if (bin_type == "square" | bin_type == "rectangle"){
        if (length(bin_param) != 2){
            stop("Invalid input bin_param, bin_param should be a vector of
                length two for rectangle/square bins")
        }
        bin_length <- bin_param[1] * bin_param[2]
    }else{
        stop("Input bin_type is not supported. Supported bin_type is
            rectangle/square or hexagon.")
    }
    #st_time = Sys.time()
    return (bin_length)
}

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
#' @param test_genes A vector of strings giving the name of the genes you want
#' to create gene vector.
#' @return outputs a list object with the following components
#' \item{trans_lst }{ A list of named dataframes. Each dataframe refers to 
#' one sample and shows the transcript detection coordinates for each gene. 
#' The name matches the input sample_names}
#' \item{cm_lst }{A list of named dataframes containing the count matrix for
#' each sample. The name matches the input sample_names} 
#' 
#' @importFrom methods is
#' @importFrom BumpyMatrix unsplitAsDataFrame
#' @importFrom SpatialExperiment molecules
convert_data <- function(x, sample_names, test_genes){
    n_samples <- length(sample_names)
    cm_lst <- vector("list", n_samples)
    names(cm_lst) <- sample_names
    trans_lst <- vector("list", n_samples)
    names(trans_lst) <- sample_names
    
    if ((is.vector(test_genes) == FALSE)){
        stop("Invalid input test_genes, should be a vector of character")}
    primary_class <- class(x)[1]
    if (primary_class == "SingleCellExperiment"){
        if (any(!sample_names %in% names(x@assays@data@listData))){
    stop("Input sample_names does not match the sample information in x.")}
        if(length(setdiff(test_genes,row.names(x))) >0){
            stop("Invalid input test_genes, can not match test_genes from x")}
        cm_lst <- x@assays@data@listData
        trans_lst <- NULL
    }else if (primary_class == "SpatialExperiment" | 
            primary_class ==  "SpatialFeatureExperiment"){
        if (!all(sample_names %in% unique(x$sample_id))){
            stop("Input sample_names does not match the 
                sample information in x.")}
        if(length(setdiff(test_genes,row.names(x))) >0){
            stop("Invalid input test_genes, can not match test_genes from x")}
        #all_cm <- x@assays@data$counts
        for (sp in sample_names){
            # x$imgData <- NULL
            sub_x <- x[,x$sample_id==sp]
            cm_lst[[sp]] <- sub_x@assays@data$counts 
        if (!is.null(x@assays@data$molecules)){
            for (sp in sample_names){
                sub_x <- x[,x$sample_id==sp]
                #transcript_mlist <- sub_x@assays@data$molecules
                tmp_data <- BumpyMatrix::unsplitAsDataFrame(molecules(sub_x))
                tmp_data <- as.data.frame(tmp_data)
                colnames(tmp_data) <- c("feature_name","cell_id","x","y")
                tmp_data$x <- as.numeric(tmp_data$x)
                tmp_data$y <- as.numeric(tmp_data$y)
#       genes_with_coords <- unique(tmp_data$feature_name)
#     if(length(intersect(test_genes,genes_with_coords)) == 0){
# stop("Invalid input test_genes, can not match test_genes from x")}
                trans_lst[[sp]]<- tmp_data
                }
        }else{ trans_lst <- NULL }
            }
    }else{
        stop(sprintf("The input class of 'x' is not supported. 
Please convert 'x' to one of the following supported types: 
SingleCellExperiment, SpatialExperiment, or SpatialFeatureExperiment.")) }
    return (list(trans_lst = trans_lst, cm_lst = cm_lst))
}

#' Vectorise the spatial coordinates
#'
#' @description
#' This function will convert the coordinates into a numeric vector for
#' genes and clusters.
#'
#' @details
#' This function can be used to generate input for \code{\link{lasso_markers}}
#' by specifying all the parameters.
#'
#' Suppose the input data contains \eqn{n} genes, \eqn{c} clusters, and
#' \eqn{k} samples, we
#' want to use \eqn{a \times a} square bin to convert the coordinates
#' of genes and clusters into 1d vectors.
#'
#' If \eqn{k=1}, the returned list will contain one matrix for gene vectors
#' (\code{gene_mt}) of dimension \eqn{a^2 \times n} and one matrix for
#' cluster vectors (\code{cluster_mt}) of dimension \eqn{a^2 \times c}.
#'
#' If \eqn{k>1}, gene and cluster vectors are constructed for each sample
#' separately and concat together. There will be additional k columns on the
#' returned \code{cluster_mt}, which is the one-hot encoding of the
#' sample information.
#'
#' Moreover, this function can vectorise genes and clusters separately based
#' on the input. If \code{x} is NULL, this function will
#' return vectorised clusters based on \code{cluster_info}.
#' If \code{cluster_info} is NULL, this function will return vectorised genes
#' based on \code{x}.
#'
#' @param x a SingleCellExperiment or SpatialExperiment or 
#' SpatialFeatureExperiment object
#' @param cluster_info A dataframe/matrix containing the centroid coordinates,
#' cluster label and sample for each cell.The column names must include
#' "x" (x coordinate), "y" (y coordinate),
#' "cluster" (cluster label) and "sample" (sample).
#' @param sample_names a vector of strings giving the sample names 
#' @param bin_type A string indicating which bin shape is to be used for
#' vectorization. One of "square" (default), "rectangle", or "hexagon".
#' @param bin_param A numeric vector indicating the size of the bin. If the
#' \code{bin_type} is "square" or "rectangle", this will be a vector of length
#' two giving the numbers of rectangular quadrats in the x and y directions. If
#' the \code{bin_type} is "hexagonal", this will be a number giving the side
#' length of hexagons. Positive numbers only.
#' @param test_genes A vector of strings giving the name of the genes you want
#' to create gene vector. This will be used as column names for one of the 
#' result matrix \code{gene_mt}.
#' @param w_x A numeric vector of length two specifying the x coordinate
#' limits of enclosing box.
#' @param w_y A numeric vector of length two specifying the y coordinate
#' limits of enclosing box.
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

#' @importFrom spatstat.geom ppp
#' @importFrom spatstat.geom quadratcount
#' @importFrom dplyr distinct
#' @importFrom stats model.matrix
#' @importFrom magrittr "%>%"
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel register
#' @importFrom BiocParallel SnowParam
#' @export
#' @examples
#' library(SpatialExperiment)
#' set.seed(100)
#' #  simulate coordinates for clusters
#' df_clA = data.frame(x = rnorm(n=100, mean=20, sd=5),
#'                  y = rnorm(n=100, mean=20, sd=5), cluster="A")
#' df_clB = data.frame(x = rnorm(n=100, mean=100, sd=5),
#'                 y = rnorm(n=100, mean=100, sd=5), cluster="B")
#'
#' clusters = rbind(df_clA, df_clB)
#' clusters$sample="sample1"
#'
#' # simulate coordinates for genes
#' trans_info = data.frame(rbind(cbind(x = rnorm(n=10, mean=20,sd=5),
#'                                 y = rnorm(n=10, mean=20, sd=5),
#'                                  feature_name="gene_A1"),
#'                            cbind(x = rnorm(n=10, mean=20, sd=5),
#'                                  y = rnorm(n=10, mean=20, sd=5),
#'                                  feature_name="gene_A2"),
#'                            cbind(x = rnorm(n=10, mean=100, sd=5),
#'                                  y = rnorm(n=10, mean=100, sd=5),
#'                                  feature_name="gene_B1"),
#'                            cbind(x = rnorm(n=10, mean=100, sd=5),
#'                                  y = rnorm(n=10, mean=100, sd=5),
#'                                  feature_name="gene_B2")))
#' trans_info$x=as.numeric(trans_info$x)
#' trans_info$y=as.numeric(trans_info$y)
#' trans_info$cell = sample(c("cell1","cell2","cell2"),replace=TRUE,
#'                         size=nrow(trans_info))
#' trans_mol <- BumpyMatrix::splitAsBumpyMatrix(
#'     trans_info[, c("x", "y")], 
#'     row = trans_info$feature_name, col = trans_info$cell )
#' spe<- SpatialExperiment(
#'      assays = list(molecules = trans_mol),sample_id ="sample1" )
#' w_x =  c(min(floor(min(trans_info$x)),
#'          floor(min(clusters$x))),
#'       max(ceiling(max(trans_info$x)),
#'           ceiling(max(clusters$x))))
#' w_y =  c(min(floor(min(trans_info$y)),
#'           floor(min(clusters$y))),
#'       max(ceiling(max(trans_info$y)),
#'           ceiling(max(clusters$y))))
#' vecs_lst = get_vectors(x=spe,sample_names=c("sample1"),
#'                     cluster_info = clusters,
#'                     bin_type = "square",
#'                     bin_param = c(5,5),
#'                     test_genes =c("gene_A1","gene_A2","gene_B1","gene_B2"),
#'                     w_x = w_x, w_y=w_y)
#'
get_vectors<- function(x, cluster_info, sample_names, bin_type, bin_param,
                        test_genes, w_x, w_y, use_cm = FALSE, n_cores=1){
    # check input
    if ((is.null(x) ==TRUE) &  (is.null(cluster_info) == TRUE) ){
        stop("Invalid input, no coordinates information is specified") }
    # convert SingleCellExperiment/SpatialExperiment/SpatialFeatureExperiment
    if (is.null(x) ==FALSE){
        re_lst <- convert_data(x, sample_names, test_genes)
        trans_lst <- re_lst$trans_lst
        cm_lst <- re_lst$cm_lst
        # decide how to define gene vectors 
        if ((is.null(trans_lst) ==FALSE) & (is.null(cm_lst) == FALSE)){
            if (use_cm==TRUE){
                trans_lst <- NULL
            }else{
                cm_lst <- NULL }}
    }else{
        trans_lst <- NULL
        cm_lst <- NULL
    }
    # will create spatial vectors for every cluster
    if ((is.null(cluster_info) == FALSE)){
        req_cols <- c("x","y","cluster","sample")
        if (length(setdiff(req_cols, colnames(cluster_info))) != 0){
            stop("Invalid columns in input clusters. Input clusters must
            contain columns 'x', 'y', 'cluster', 'sample'
                for every cell")}
        if ((is.null(cm_lst) == FALSE)){
            if (length(setdiff(unique(cluster_info$sample),names(cm_lst)))!=0){
                stop("Mismatched sample names in cluster_info and x") }
            # must contain cell id
            req_cols <- c("x","y","cluster","sample","cell_id")
            if (length(setdiff(req_cols, colnames(cluster_info))) !=0){
                stop("Invalid columns in input clusters. Input cluster_info
                    must contain columns 'x', 'y', 'cluster',
                    'sample','cell_id' for every cell") } } }
    # must provide cluster info if gene vectors are created from count matrix
    if ((is.null(cm_lst) == FALSE) & (is.null(cluster_info) == TRUE)){
        stop("Missing cluster information to build gene vector matrix.")}
    bin_length<-check_binning(bin_param=bin_param,bin_type=bin_type,
                        w_x=w_x,w_y=w_y)
    # register for BiocParallel 
    register(SnowParam(workers = n_cores, type = "SOCK"))
    # with cluster information
    if (is.null(cluster_info) == FALSE){
        vec_cluster <- get_cluster_vectors(cluster_info=cluster_info,
                                    bin_length=bin_length, bin_type=bin_type,
                                    bin_param=bin_param,w_x=w_x,w_y=w_y,
                                    sample_names=sample_names) }
    # with gene information
    if (is.null(trans_lst) == FALSE & is.null(cm_lst)== TRUE){
        vec_gene_mt<-get_gene_vectors_tr(trans_lst=trans_lst, 
                                test_genes=test_genes,bin_type=bin_type, 
                                bin_param=bin_param,
                                bin_length=bin_length,w_x=w_x,w_y=w_y)}
    if (is.null(trans_lst) == TRUE & is.null(cm_lst)==FALSE &
        is.null(cluster_info)==FALSE){
        # use count matrix to build gene vector matrix
        vec_gene_mt<-get_gene_vectors_cm(cluster_info=cluster_info,
                        cm_lst=cm_lst, bin_type=bin_type, bin_param=bin_param,
                        test_genes=test_genes,w_x=w_x, w_y=w_y) }
    result <- list()
    if ((is.null(cluster_info) == FALSE) ){
        result$cluster_mt <- as.matrix(vec_cluster)}
    if ((is.null(trans_lst) == FALSE) | (is.null(cm_lst) == FALSE)){
        result$gene_mt <- as.matrix(vec_gene_mt)}
    return (result)
}
