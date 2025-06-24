
#' Create a marker gene result object for correlation approach
#'
#' This function creates a structured output object named 'cor_mg_result' for 
#' storing the permutation results. 
#' The object contains three matrices:
#' @param obs.stat A matrix containing the correlation coefficients for 
#' each pair of genes and cluster vectors.
#' @param perm.pval A matrix containing the raw permutation p-value for 
#' each pair of genes and cluster.
#' @param perm.pval.adj A matrix containing the adjusted permutation p-value
#' for each pair of genes and cluster.
#' @return An S3 object of class 'cor_mg_result' which includes three 
#' matrices.

.create_cor_mg_result <- function(obs.stat, perm.pval, perm.pval.adj) {
    structure(list(obs.stat = obs.stat, perm.pval = perm.pval, 
                    perm.pval.adj=perm.pval.adj), class = "cor_mg_result")
}

#' Get observed correlation cor_mg_result
#'
#' Accessor function to retrieve the observed correlation from 
#' an 'cor_mg_result' object.
#' @param obj An 'cor_mg_result' object.
#' @return A matrix contains the observation statistic for
#' every gene and every cluster. Each row refers to a gene, and each column
#' refers to a cluster
#' @export
#' @examples 
#' library(SpatialExperiment)
#' library(BumpyMatrix)
#' set.seed(100)
#' # simulate coordinates for clusters
#' df_clA <- data.frame(x = rnorm(n=10, mean=20, sd=5),
#'                     y = rnorm(n=10, mean=20, sd=5), cluster="A")
#' df_clB <- data.frame(x = rnorm(n=10, mean=100, sd=5),
#'                     y = rnorm(n=10, mean=100, sd=5), cluster="B")
#' clusters <- rbind(df_clA, df_clB)
#' clusters$sample="sample1"
#' # simulate coordinates for genes
#' trans_info <- data.frame(rbind(cbind(x = rnorm(n=10, mean=20, sd=5),
#'                                     y = rnorm(n=10, mean=20, sd=5),
#'                                     feature_name="gene_A1"),
#'                     cbind(x = rnorm(n=10, mean=20, sd=5),
#'                                     y = rnorm(n=10, mean=20, sd=5),
#'                                     feature_name="gene_A2"),
#'                     cbind(x = rnorm(n=10, mean=100, sd=5),
#'                                     y = rnorm(n=10, mean=100, sd=5),
#'                                     feature_name="gene_B1"),
#'                     cbind(x = rnorm(n=10, mean=100, sd=5),
#'                                     y = rnorm(n=10, mean=100, sd=5),
#'                                     feature_name="gene_B2")))
#' trans_info$x<-as.numeric(trans_info$x)
#' trans_info$y<-as.numeric(trans_info$y)
#' trans_info$cell =  rep(paste("cell",1:20, sep=""), times=2)
#' mol <- BumpyMatrix::splitAsBumpyMatrix(
#'      trans_info[, c("x", "y")], 
#'      row = trans_info$feature_name, col = trans_info$cell )
#' spe_sample1 <- SpatialExperiment(
#'         assays = list(molecules = mol),sample_id ="sample1" )
#' set.seed(100)
#' corr_res <- compute_permp(x=spe_sample1,
#'              cluster_info=clusters,
#'              perm.size=10,
#'              bin_type="square",
#'              bin_param=c(2,2),
#'              test_genes=unique(trans_info$feature_name),
#'              correlation_method = "pearson",
#'              n_cores=1,
#'              correction_method="BH")
#' # observed correlation 
#' obs_corr <- get_cor(corr_res)
get_cor <- function(obj) {
    stopifnot(inherits(obj, "cor_mg_result"))
    return(obj$obs.stat)
}


#' Get permutation p value from cor_mg_result
#'
#' Accessor function to retrieve the raw permutation p-value from 
#' an 'cor_mg_result' object.
#' @param obj An 'cor_mg_result' object.
#' @return A matrix contains the raw permutation p-value.
#' Each row refers to a gene, and each column refers to a cluster.
#' @export
#' @examples 
#' library(SpatialExperiment)
#' library(BumpyMatrix)
#' set.seed(100)
#' # simulate coordinates for clusters
#' df_clA <- data.frame(x = rnorm(n=10, mean=20, sd=5),
#'                     y = rnorm(n=10, mean=20, sd=5), cluster="A")
#' df_clB <- data.frame(x = rnorm(n=10, mean=100, sd=5),
#'                     y = rnorm(n=10, mean=100, sd=5), cluster="B")
#' clusters <- rbind(df_clA, df_clB)
#' clusters$sample<-"sample1"
#' # simulate coordinates for genes
#' trans_info <- data.frame(rbind(cbind(x = rnorm(n=10, mean=20, sd=5),
#'                                     y = rnorm(n=10, mean=20, sd=5),
#'                                     feature_name="gene_A1"),
#'                     cbind(x = rnorm(n=10, mean=20, sd=5),
#'                                     y = rnorm(n=10, mean=20, sd=5),
#'                                     feature_name="gene_A2"),
#'                     cbind(x = rnorm(n=10, mean=100, sd=5),
#'                                     y = rnorm(n=10, mean=100, sd=5),
#'                                     feature_name="gene_B1"),
#'                     cbind(x = rnorm(n=10, mean=100, sd=5),
#'                                     y = rnorm(n=10, mean=100, sd=5),
#'                                     feature_name="gene_B2")))
#' trans_info$x<-as.numeric(trans_info$x)
#' trans_info$y<-as.numeric(trans_info$y)
#' trans_info$cell =  rep(paste("cell",1:20, sep=""), times=2)
#' mol <- BumpyMatrix::splitAsBumpyMatrix(
#'      trans_info[, c("x", "y")], 
#'      row = trans_info$feature_name, col = trans_info$cell )
#' spe_sample1 <- SpatialExperiment(
#'         assays = list(molecules = mol),sample_id ="sample1" )

#' set.seed(100)
#' corr_res <- compute_permp(x=spe_sample1,
#'              cluster_info=clusters,
#'              perm.size=10,
#'              bin_type="square",
#'              bin_param=c(2,2),
#'              test_genes=unique(trans_info$feature_name),
#'              correlation_method = "pearson",
#'              n_cores=1,
#'              correction_method="BH")
#'              
#' # raw permutation p-value
#' perm_p <- get_perm_p(corr_res)
get_perm_p <- function(obj) {
    stopifnot(inherits(obj, "cor_mg_result"))
    return(obj$perm.pval)
}

#' Get permutation adjusted p value from cor_mg_result
#'
#' Accessor function to retrieve the permutation adjusted p-value from 
#' an 'cor_mg_result' object.
#' @param obj An 'cor_mg_result' object.
#' @return A matrix contains the adjusted permutation
#' p-value. Each row refers to a gene, and each column refers to a cluster.
#' @export
#' @examples 
#' library(SpatialExperiment)
#' library(BumpyMatrix)
#' set.seed(100)
#' # simulate coordinates for clusters
#' df_clA <- data.frame(x = rnorm(n=10, mean=20, sd=5),
#'                     y = rnorm(n=10, mean=20, sd=5), cluster="A")
#' df_clB <- data.frame(x = rnorm(n=10, mean=100, sd=5),
#'                     y = rnorm(n=10, mean=100, sd=5), cluster="B")
#' clusters <- rbind(df_clA, df_clB)
#' clusters$sample="sample1"
#' # simulate coordinates for genes
#' trans_info <- data.frame(rbind(cbind(x = rnorm(n=10, mean=20, sd=5),
#'                                     y = rnorm(n=10, mean=20, sd=5),
#'                                     feature_name="gene_A1"),
#'                     cbind(x = rnorm(n=10, mean=20, sd=5),
#'                                     y = rnorm(n=10, mean=20, sd=5),
#'                                     feature_name="gene_A2"),
#'                     cbind(x = rnorm(n=10, mean=100, sd=5),
#'                                     y = rnorm(n=10, mean=100, sd=5),
#'                                     feature_name="gene_B1"),
#'                     cbind(x = rnorm(n=10, mean=100, sd=5),
#'                                     y = rnorm(n=10, mean=100, sd=5),
#'                                     feature_name="gene_B2")))
#' trans_info$x<-as.numeric(trans_info$x)
#' trans_info$y<-as.numeric(trans_info$y)
#' trans_info$cell =  rep(paste("cell",1:20, sep=""), times=2)
#' mol <- BumpyMatrix::splitAsBumpyMatrix(
#'      trans_info[, c("x", "y")], 
#'      row = trans_info$feature_name, col = trans_info$cell )
#' spe_sample1 <- SpatialExperiment(
#'         assays = list(molecules = mol),sample_id ="sample1" )

#' set.seed(100)
#' corr_res <- compute_permp(x=spe_sample1,
#'              cluster_info=clusters,
#'              perm.size=10,
#'              bin_type="square",
#'              bin_param=c(2,2),
#'              test_genes=unique(trans_info$feature_name),
#'              correlation_method = "pearson",
#'              n_cores=1,
#'              correction_method="BH")
#' # adjusted permutation p-value
#' adjusted_perm_p <- get_perm_adjp(corr_res)
get_perm_adjp <- function(obj) {
    stopifnot(inherits(obj, "cor_mg_result"))
    return(obj$perm.pval.adj)
}