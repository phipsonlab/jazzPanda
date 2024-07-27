#' helper function to check the inputs passed to create geneset function
#' 
#' @param data_lst A list of named matrices containing the coordinates of 
#' transcripts.
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
#' @param cluster_info A dataframe/matrix containing the centroid coordinates,
#' cluster sample label for each cell.The column names must include
#' "x" (x coordinate), "y" (y coordinate),
#' "cluster" (cluster label) and "sample" (sample).
#' @return A list of two elements

check_geneset_input <- function(data_lst, bin_type, bin_param, 
                                    w_x, w_y, cluster_info){
    # binning
    bin_length <- 0
    if (bin_type == "hexagon"){
        if (length(bin_param) != 1){
            stop("Invalid input bin_param, bin_param should be a vector of 
                length 2 for hexagon bins")
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
    
    if (length(unique(vapply(data_lst,class,FUN.VALUE = character(1)))) != 1) {
        stop("The input data_lst contains elements from multiple classes")
    }
    
    if (length(names(data_lst)) == 0){
        stop("data_lst should be a named list")
    }
    
    use_cm <- FALSE
    for (rpp in data_lst){
        if (typeof(rpp) == "list"){
            req_cols <- c('x', 'y', 'feature_name' )
            if (length(setdiff(req_cols, colnames(rpp$trans_info))) > 0){
                stop("Invalid column names detected in input data_lst. 
                Must contain columns 'x', 'y', 
                'feature_name' for every transcript")}
        }else if (is.matrix(rpp) ==TRUE | is.data.frame(rpp)==TRUE){
            if (is.null(cluster_info)==TRUE){
                stop("Missing cluster information to build gene vector matrix")
            }
            # must contain cell id
            req_cols <- c("x","y","cluster","sample","cell_id")
            if (length(setdiff(req_cols, colnames(cluster_info))) !=0){
                stop("Invalid columns in input clusters. Input cluster_info 
            must contain columns 'x', 'y', 'cluster', 
            'sample','cell_id' for every cell")
            }
            use_cm <- TRUE
        }
        
    }
    return (list(bin_length = bin_length, use_cm=use_cm))
}

#' Convert the coordinates of set of genes into vectors.
#'
#' @param data_lst A list of named matrices containing the coordinates of 
#' transcripts.
#' @param name_lst A named list of strings giving the name of features that are
#' treated as background.
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
#' @param cluster_info A dataframe/matrix containing the centroid coordinates,
#' cluster and sample label for each cell.The column names must include
#' "x" (x coordinate), "y" (y coordinate),
#' "cluster" (cluster label) and "sample" (sample).
#'
#' @importFrom stats setNames
#' @importFrom magrittr "%>%"
#' @return A matrix contains the sum count in each grid.
#' Each row refers to a grid, each column refers to a set in \code{name_lst}.
#' The column name will match the names in \code{name_lst}.
#' @export
#'
#' @examples
#' set.seed(15)
#' trans = as.data.frame(rbind(cbind(x = runif(10, min=1, max=10),
#'                                 y = runif(10, min=1, max=10),
#'                                 feature_name="A"),
#'                          cbind(x = runif(5, min=10, max=24),
#'                                y = runif(5, min=1, max=10),
#'                                feature_name="B"),
#'                          cbind(x = runif(10, min=10, max=24),
#'                                y = runif(10, min=10, max=24),
#'                                feature_name="C")))
#' trans$x = as.numeric(trans$x)
#' trans$y = as.numeric(trans$y)
#' data=list(trans_info=trans)
#' geneset_res = create_genesets(data_lst=list("rep1"= data),
#'                            name_lst=list(dummy_A=c("A","C"),
#'                                          dummy_B=c("A","B","C")),
#'                            bin_type="square",
#'                            bin_param=c(2,2),
#'                            w_x=c(0,25), w_y=c(0,25))
#'
create_genesets<-function(data_lst, name_lst, bin_type="square",
                            bin_param, w_x=w_x, w_y=w_y,
                            cluster_info){
    
    input_info <- check_geneset_input(data_lst=data_lst, bin_type=bin_type, 
                        bin_param=bin_param, w_x=w_x, w_y=w_y, 
                        cluster_info=cluster_info)
    use_cm <- input_info$use_cm
    bin_length <- input_info$bin_length
    if (bin_type == "hexagon"){
        w <- owin(xrange=w_x, yrange=w_y)
        H <- hextess(W=w, bin_param[1])
    }
    n_samples <- length(data_lst)
    vec_background <- as.data.frame(matrix(0, ncol=length(name_lst),
                            nrow=bin_length*n_samples))
    colnames(vec_background) <- names(name_lst)

    # iterate over each set over the name_lst
    for (nm in names(name_lst)){
        vec_name<- c()
        for (rp_nm in names(data_lst)){
            rpp <- data_lst[[rp_nm]]
            if (use_cm==TRUE){
                cm_lst <- setNames(list(rpp), rp_nm)
                # use cell-level coordinates and count matrix
                vec_g <- get_gene_vectors_cm(cluster_info=cluster_info,
                                            cm_lst=cm_lst, bin_type=bin_type,
                                            bin_param=bin_param,
                                            all_genes=name_lst[[nm]],
                                            w_x=w_x, w_y=w_y)
                vec_name <- c(vec_name, rowSums(vec_g))
            }else{
                curr <- rpp$trans_info[rpp$trans_info$feature_name %in% 
                                name_lst[[nm]],
                                c("x", "y")] %>% distinct()
                gene_ppp <- ppp(curr$x,curr$y,w_x, w_y)

                if (bin_type == "hexagon"){
                    vec_g <- as.vector(t(quadratcount(gene_ppp, tess=H)))
                    vec_name <- c(vec_name, vec_g)
                }else{
                    vec_g <- as.vector(t(quadratcount(gene_ppp, bin_param[1],
                                                        bin_param[2])))
                    vec_name<- c(vec_name, vec_g)
                }
            }
        }
        vec_background[,nm] <- vec_name
    }
    return (vec_background)
}
