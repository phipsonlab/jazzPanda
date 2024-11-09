#' Compute observation statistic for permutation framework
#'
#' @param data A list of matrices containing the coordinates of transcripts.
#' @param cluster_info A dataframe/matrix containing the centroid
#' coordinates and cluster label for each cell.The column names should
#' include "x" (x coordinate), "y" (y coordinate),
#' and "cluster" (cluster label).
#' @param correlation_method A parameter pass to \code{\link{cor}},
#' indicating which correlation coefficient is to be computed.
#' One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param all_genes A vector of strings giving the name of the genes you want to
#' test correlation for.
#' @param bin_type A string indicating which bin shape is to be used for
#' vectorization. One of "square" (default), "rectangle", or "hexagon".
#' @param bin_param A numeric vector indicating the size of the bin. If the
#' \code{bin_type} is "square" or "rectangle", this will be a vector of length
#' two giving the numbers of rectangular quadrats in the x and y directions. If
#' the \code{bin_type} is "hexagonal", this will be a number giving the side
#' length of hexagons. Positive numbers only.
#' @param w_x a numeric vector of length two specifying the x coordinate
#' limits of enclosing box.
#' @param w_y a numeric vector of length two specifying the y coordinate
#' limits of enclosing box.
#'
#' @importFrom stats cor
#' @return A named list with the following components
#' \item{\code{obs.stat}  }{ A matrix contains the observation statistic for
#' every gene and every cluster. Each row refers to a gene, and each column
#' refers to a cluster}
#' \item{\code{gene_mt}  }{ contains the transcript count in each grid.
#' Each row refers to a grid, and each column refers to a gene.}
compute_observation<- function(data, cluster_info,
                                correlation_method = "pearson",
                                all_genes=all_genes,
                                bin_type, bin_param, w_x, w_y){
    
    vectors_lst <- get_vectors(trans_lst=data,
                                cluster_info = cluster_info,
                                bin_type=bin_type,
                                bin_param=bin_param,
                                all_genes=all_genes,
                                w_x=w_x, w_y=w_y)

    # calculate correation of permuted clusters and gene
    obs.stat <- cor(x=vectors_lst$gene_mt, y=vectors_lst$cluster_mt,
                        method=correlation_method)
    return (list(obs.stat = obs.stat, gene_mt=vectors_lst$gene_mt))
}

#' Compute permutation statistics for permutation framework
#'
#' @param cluster_info A dataframe/matrix containing the centroid
#' coordinates and cluster label for each cell.The column names
#' should include "x" (x coordinate), "y" (y coordinate), and
#' "cluster" (cluster label).
#' @param perm.size A positive number specifying permutation times
#' @param correlation_method A parameter pass to \code{\link{cor}},
#' indicating which correlation coefficient is to be computed.
#' One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param bin_type A string indicating which bin shape is to be used for
#' vectorization. One of "square" (default), "rectangle", or "hexagon".
#' @param bin_param A numeric vector indicating the size of the bin. If the
#' \code{bin_type} is "square" or "rectangle", this will be a vector of length
#' two giving the numbers of rectangular quadrats in the x and y directions. If
#' the \code{bin_type} is "hexagonal", this will be a number giving the side
#' length of hexagons. Positive numbers only.
#' @param n.cores A positive number specifying number of cores used for
#' parallelizing permutation testing. Default is one core
#' (sequential processing).
#' @param w_x a numeric vector of length two specifying the x coordinate
#' limits of enclosing box.
#' @param w_y a numeric vector of length two specifying the y coordinate
#' limits of enclosing box.
#' @param gene_mt A matrix contains the transcript count in each grid.
#' Each row refers to a grid, and each column refers to a gene.
#' @param cluster_names A list of strings giving the name and order of the
#' clusters
#' @importFrom foreach foreach
#' @importFrom foreach `%dopar%`
#' @importFrom stats cor
#' @return A matrix with permutation statistics
#'
compute_permutation<- function(cluster_info, perm.size = 1000,
                                correlation_method = "pearson",  bin_type,
                                bin_param, n.cores=1, w_x,w_y, gene_mt,
                                cluster_names){
    n_clusters <- length(unique(cluster_info$cluster))
    t.perm.array<- array(0, dim = c(ncol(gene_mt),
                            length(cluster_names),perm.size))
    if (bin_type == "hexagon"){
        if (length(bin_param) != 1){
            stop("Invalid input bin_param, bin_param should be a
                    vector of length 2 for hexagon bins") }
        w <-owin(xrange=w_x, yrange=w_y)
        H <-hextess(W=w, bin_param[1])
        bin_length <- length(H$tiles)
    }else if (bin_type == "square" | bin_type == "rectangle"){
    if (length(bin_param) != 2){
        stop("Invalid input bin_param, bin_param should be
                a vector of length 2 for rectangle/square bins") }
    bin_length <- bin_param[1] * bin_param[2]
    }else{
        stop("Input bin_type is not supported.
                Supported bin_type is rectangle/square or hexagon.") }
    my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK" )
    doParallel::registerDoParallel(cl = my.cluster)
    i <- NULL
    t.perm.array<-foreach (i = seq_len(perm.size)) %dopar% {
        cell_cluster <- cluster_info
        # permutate the cluster labels
        cell_cluster$cluster <- sample(cell_cluster$cluster,
                                        size = length(cell_cluster$cluster),
                                        replace = FALSE)
        cluster_mt <- matrix(0, ncol=length(cluster_names), nrow=bin_length)
        colnames(cluster_mt) <- cluster_names
        # a matrix of cluster vector, each column as a vector for a cluster
        for (i_cluster in cluster_names){
            x_loc <- cell_cluster[cell_cluster$cluster==i_cluster, "x"]
            y_loc <- cell_cluster[cell_cluster$cluster==i_cluster, "y"]
            cluster_ppp <- spatstat.geom::ppp(x_loc,y_loc,w_x, w_y)
            if (bin_type == "hexagon"){
                cm_cluster <- spatstat.geom::quadratcount(cluster_ppp, tess=H)
            }else{
                cm_cluster <- spatstat.geom::quadratcount(cluster_ppp, 
                                                        bin_param[1],
                                                        bin_param[2]) }
            cluster_mt[,i_cluster] <- as.vector(t(cm_cluster)) }
        # calculate correlation of permuted clusters and gene
        perm.stat <- cor(x=gene_mt, y=cluster_mt, method=correlation_method)
        t.perm.array[,,i]<- perm.stat }
    perm.array<- simplify2array(t.perm.array)
    parallel::stopCluster(cl = my.cluster)
    return (list(t.perm = perm.array))
}


#' Calculate a p-value for correlation with permutation.

#' @description
#' This function will run permutation framework to compute a p-value for the
#' correlation between the vectorised genes and clusters each cluster for one 
#' sample.
#'
#' @details
#' To get a permutation p-value for the correlation between a gene
#' and a cluster, this function will permute the cluster label for
#' each cell randomly, and calculate correlation between the genes and
#' permuted clusters. This process will be repeated for \code{perm.size}
#' times, and permutation p-value is calculated as the probability of
#' permuted correlations larger than the observation correlation.
#'
#' @param data A list of matrices containing the coordinates of transcripts.
#' @param cluster_info A dataframe/matrix containing the centroid coordinates
#' and cluster label for each cell.The column names should include "x"
#' (x coordinate), "y" (y coordinate), and "cluster" (cluster label).
#' @param perm.size A positive number specifying permutation times
#' @param bin_type A string indicating which bin shape is to be used for
#' vectorization. One of "square" (default), "rectangle", or "hexagon".
#' @param bin_param A numeric vector indicating the size of the bin. If the
#' \code{bin_type} is "square" or "rectangle", this will be a vector of length
#' two giving the numbers of rectangular quadrats in the x and y directions. If
#' the \code{bin_type} is "hexagonal", this will be a number giving the side
#' length of hexagons. Positive numbers only.
#' @param all_genes A vector of strings giving the name of the genes you
#' want to test correlation for.
#' \code{gene_mt}.
#' @param correlation_method A parameter pass to \code{\link{cor}} indicating
#' which correlation coefficient is to be computed.
#' One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param n.cores A positive number specifying number of cores used for
#' parallelizing permutation testing. Default is one core
#' (sequential processing).
#' @param w_x a numeric vector of length two specifying the x coordinate
#' limits of enclosing box.
#' @param w_y a numeric vector of length two specifying the y coordinate
#' limits of enclosing box.
#'
#' @param correction_method A character string pass to \code{\link{p.adjust}}
#' specifying the correction method for multiple testing .
#'
#' @return A named list with the following components
#' \item{\code{obs.stat}  }{ A matrix contains the observation statistic for
#' every gene and every cluster. Each row refers to a gene, and each column
#' refers to a cluster}
#' \item{\code{perm.arrays}  }{ A three dimensional array.
#' The first two dimensions represent the correlation between the genes and
#' permuted clusters. The third dimension refers to the different permutation
#' runs. }
#' \item{\code{perm.pval}  }{A matrix contains the raw permutation p-value.
#' Each row refers to a gene, and each column refers to a cluster}
#' \item{\code{perm.pval.adj}  }{A matrix contains the adjusted permutation
#' p-value. Each row refers to a gene, and each column refers to a cluster}

#' @importFrom stats p.adjust
#' @importFrom foreach foreach
#' @importFrom foreach `%dopar%`
#' @export
#' @examples
#'
#' set.seed(100)
#' # simulate coordinates for clusters
#' df_clA = data.frame(x = rnorm(n=100, mean=20, sd=5),
#'                    y = rnorm(n=100, mean=20, sd=5), cluster="A")
#' df_clB = data.frame(x = rnorm(n=100, mean=100, sd=5),
#'                   y = rnorm(n=100, mean=100, sd=5), cluster="B")

#' clusters = rbind(df_clA, df_clB)
#' clusters$sample="rep1"
#' # simulate coordinates for genes
#' trans_info = data.frame(rbind(cbind(x = rnorm(n=100, mean=20, sd=5),
#'                                     y = rnorm(n=100, mean=20, sd=5),
#'                                  feature_name="gene_A1"),
#'                            cbind(x = rnorm(n=100, mean=20, sd=5),
#'                                  y = rnorm(n=100, mean=20, sd=5),
#'                                  feature_name="gene_A2"),
#'                            cbind(x = rnorm(n=100, mean=100, sd=5),
#'                                  y = rnorm(n=100, mean=100, sd=5),
#'                                  feature_name="gene_B1"),
#'                            cbind(x = rnorm(n=100, mean=100, sd=5),
#'                                  y = rnorm(n=100, mean=100, sd=5),
#'                                  feature_name="gene_B2")))
#' trans_info$x=as.numeric(trans_info$x)
#' trans_info$y=as.numeric(trans_info$y)

#' w_x =  c(min(floor(min(trans_info$x)),
#'             floor(min(clusters$x))),
#'         max(ceiling(max(trans_info$x)),
#'             ceiling(max(clusters$x))))
#' w_y =  c(min(floor(min(trans_info$y)),
#'          floor(min(clusters$y))),
#'       max(ceiling(max(trans_info$y)),
#'           ceiling(max(clusters$y))))
#'           
#' perm_res_lst = compute_permp(data=list("sample1"=trans_info),
#'                     cluster_info=clusters,
#'                     perm.size=100,
#'                     bin_type="square",
#'                     bin_param=c(2,2),
#'                     all_genes=unique(trans_info$feature_name),
#'                     correlation_method = "pearson",
#'                     n.cores=2,
#'                     correction_method="BH",
#'                     w_x=w_x ,
#'                     w_y=w_y)
#' perm_pvalue = perm_res_lst$perm.pval.adj
compute_permp<-function(data, cluster_info, perm.size, bin_type,
                        bin_param,all_genes,
                        correlation_method = "pearson", n.cores=1,
                        correction_method="BH",w_x ,w_y){
    message(sprintf("Correlation Method = %s", correlation_method))

    tm1 <- system.time(
    {
        obs_res<- compute_observation(data=data, cluster_info=cluster_info,
                                        correlation_method = correlation_method,
                                        bin_type=bin_type,all_genes=all_genes,
                                        bin_param=bin_param, w_x=w_x, w_y = w_y)
    })

    obs.stat<- obs_res$obs.stat
    if (n.cores>1 ){
        message(sprintf("Running %s permutation with %s cores in parallel",
                        perm.size, n.cores))
    }else{
        message(sprintf("Running %s permutation in sequential", perm.size))
    }


    # permutation stats
    perm_stat <- compute_permutation(cluster_info= cluster_info,
                                        perm.size = perm.size,
                                        correlation_method = correlation_method,
                                        bin_type=bin_type,
                                        bin_param=bin_param, n.cores=n.cores,
                                        w_x=w_x, w_y = w_y,
                                        gene_mt = obs_res$gene_mt,
                                        cluster_names =colnames(obs.stat) )
    # permutation stats
    perm.arrays<- perm_stat$t.perm

    perm.pvals<-apply(expand.grid(x = seq_len(dim(perm.arrays)[1]),
                                    y = seq_len(dim(perm.arrays)[2])), 1,
function(r) (sum(perm.arrays[r[1],r[2],]>obs.stat[r[1],r[2]])+1)/(perm.size+1))


    perm.pval<- matrix(perm.pvals,
                        nrow=dim(perm.arrays)[1],
                        ncol = dim(perm.arrays)[2],
                        byrow=FALSE)
    rownames(perm.pval) <- row.names(obs.stat)
    colnames(perm.pval) <- colnames(obs.stat)

    # multiple testing adjustment
    perm.pval.adj<- apply(perm.pval, 2, p.adjust, method = correction_method)
    perm.pval.adj<- as.data.frame(perm.pval.adj)
    rownames(perm.pval.adj) <- row.names(obs.stat)
    colnames(perm.pval.adj) <- colnames(obs.stat)

    return(list(obs.stat= obs.stat,
                perm.arrays=perm.arrays,
                perm.pval = perm.pval,
                perm.pval.adj=perm.pval.adj
    ))
}
