#' Create a marker gene result object for linear modelling approach
#'
#' This function creates a structured output object named 'glm_mg_result' for 
#' storing the marker gene results. 
#' The object contains two data frames: top results and full results.
#'
#' @param top_result A data frame containing top results.
#' @param full_result A data frame containing full results.
#' @return An S3 object of class 'glm_mg_result' which includes both 
#' results data frames.
.create_lm_mg_result <- function(top_result, full_result) {

    full_cols <- c("gene","cluster","glm_coef","p_value","pearson",
                    "max_gg_corr","max_gc_corr")
    top_cols <- c("gene","top_cluster","glm_coef","pearson",
                    "max_gg_corr","max_gc_corr")
    if (!inherits(top_result, "data.frame") || 
        !inherits(full_result, "data.frame")) {
        stop("Both top_result and full_result must be data frames.")
    }
    if (length(setdiff(top_cols,colnames(top_result)))>0) {
        stop("The top marker gene result is missing required columns.")
    }
    if (length(setdiff(full_cols,colnames(full_result)))>0) {
        stop("The full marker gene result is missing required columns.")
    }
    structure(list(top_result = top_result, full_result = full_result),
                    class = "glm_mg_result")
}


#' Get top lasso result from glm_mg_result
#'
#' Accessor function to retrieve the 'top_result' dataframe from 
#' an 'glm_mg_result' object.
#'
#' @param obj An 'glm_mg_result' object.
#' @param coef_cutoff A positive number giving the coefficient cutoff value.
#' Genes whose top cluster showing a coefficient value smaller than the cutoff
#' will be marked as non-marker genes ("NoSig"). Default is 0.05.
#' @return A data frame with detailed information for
#' each gene and the most relevant cluster label.
#' \itemize{
#' \item{\code{gene} Gene name}
#' \item{\code{top_cluster} The name of the most relevant cluster
#' after thresholding the coefficients. }
#' \item{\code{glm_coef} The coefficient of the selected cluster in the
#' generalised linear model.}
#' \item{\code{pearson} Pearson correlation between the gene vector and the
#' selected cluster vector. }
#' \item{\code{max_gg_corr} A number showing the maximum pearson correlation
#' for this gene vector and all other gene vectors in the input \code{gene_mt}}
#' \item{\code{max_gc_corr} A number showing the maximum pearson correlation
#' for this gene vector and every cluster vectors in the input
#' \code{cluster_mt}}
#' }
#' @export
#' @examples 
#' library(SpatialExperiment)
#' set.seed(100)
#' #  simulate coordinates for clusters
#' df_clA <- data.frame(x = rnorm(n=100, mean=20, sd=5),
#'                  y = rnorm(n=100, mean=20, sd=5), cluster="A")
#' df_clB <- data.frame(x = rnorm(n=100, mean=100, sd=5),
#'                 y = rnorm(n=100, mean=100, sd=5), cluster="B")
#'
#' clusters <- rbind(df_clA, df_clB)
#' clusters$sample<-"sample1"
#'
#' # simulate coordinates for genes
#' trans_info <- data.frame(rbind(cbind(x = rnorm(n=100, mean=20,sd=5),
#'                                 y = rnorm(n=100, mean=20, sd=5),
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
#' trans_info$x<-as.numeric(trans_info$x)
#' trans_info$y<-as.numeric(trans_info$y)
#' trans_info$cell<-sample(c("cell1","cell2","cell2"),replace=TRUE,
#'                         size=nrow(trans_info))
#' trans_mol <- BumpyMatrix::splitAsBumpyMatrix(
#'     trans_info[, c("x", "y")], 
#'     row = trans_info$feature_name, col = trans_info$cell )
#' spe<- SpatialExperiment(
#'      assays = list(molecules = trans_mol),sample_id ="sample1" )
#' w_x <- c(min(floor(min(trans_info$x)),
#'          floor(min(clusters$x))),
#'       max(ceiling(max(trans_info$x)),
#'           ceiling(max(clusters$x))))
#' w_y <- c(min(floor(min(trans_info$y)),
#'           floor(min(clusters$y))),
#'       max(ceiling(max(trans_info$y)),
#'           ceiling(max(clusters$y))))
#' vecs_lst <- get_vectors(x=spe,sample_names=c("sample1"),
#'                     cluster_info = clusters,
#'                     bin_type = "square",
#'                     bin_param = c(20,20),
#'                     test_genes =c("gene_A1","gene_A2","gene_B1","gene_B2"),
#'                     w_x = w_x, w_y=w_y)
#' lasso_res <- lasso_markers(gene_mt=vecs_lst$gene_mt,
#'                         cluster_mt = vecs_lst$cluster_mt,
#'                         sample_names=c("sample1"),
#'                         keep_positive=TRUE,
#'                         background=NULL)
#' # the top result
#' top_result<- get_top_mg(lasso_res, coef_cutoff=0.05)
get_top_mg <- function(obj, coef_cutoff=0.05) {
    stopifnot(inherits(obj, "glm_mg_result"))
    t_res <-  obj$top_result
    t_res[abs(t_res$glm_coef)<= coef_cutoff,"top_cluster"]<-"NoSig"
    t_res[abs(t_res$glm_coef)<= coef_cutoff,"glm_coef"] <- 0
    t_res[abs(t_res$glm_coef)<= coef_cutoff,"pearson"] <- 0
    return(t_res)
}

#' Get full lasso result from glm_mg_result
#'
#' Accessor function to retrieve the 'full_result' dataframe from 
#' an 'glm_mg_result' object.
#' @param obj An 'glm_mg_result' object.
#' @param coef_cutoff A positive number giving the coefficient cutoff value.
#' Genes whose cluster showing a coefficient value smaller than the cutoff
#' will be removed. Default is 0.
#' @return A data frame with detailed information for
#' each gene and the most relevant cluster label.
#'
#' \itemize{
#' \item{\code{gene} Gene name}
#' \item{\code{cluster} The name of the significant cluster after }
#' \item{\code{glm_coef} The coefficient of the selected cluster
#' in the generalised linear model.}
#' \item{\code{pearson} Pearson correlation between the gene vector and the
#' selected cluster vector. }
#' \item{\code{max_gg_corr} A number showing the maximum pearson correlation
#' for this gene vector and all other gene vectors in the input \code{gene_mt}}
#' \item{\code{max_gc_corr} A number showing the maximum pearson correlation
#' for this gene vector and every cluster vectors in the input
#' \code{cluster_mt}}
#' }
#' @export
#' @examples 
#' library(SpatialExperiment)
#' set.seed(100)
#' #  simulate coordinates for clusters
#' df_clA <- data.frame(x = rnorm(n=100, mean=20, sd=5),
#'                  y = rnorm(n=100, mean=20, sd=5), cluster="A")
#' df_clB <- data.frame(x = rnorm(n=100, mean=100, sd=5),
#'                 y = rnorm(n=100, mean=100, sd=5), cluster="B")
#'
#' clusters <- rbind(df_clA, df_clB)
#' clusters$sample<-"sample1"
#'
#' # simulate coordinates for genes
#' trans_info <- data.frame(rbind(cbind(x = rnorm(n=100, mean=20,sd=5),
#'                                 y = rnorm(n=100, mean=20, sd=5),
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
#' trans_info$x<-as.numeric(trans_info$x)
#' trans_info$y<-as.numeric(trans_info$y)
#' trans_info$cell<-sample(c("cell1","cell2","cell2"),replace=TRUE,
#'                         size=nrow(trans_info))
#' trans_mol <- BumpyMatrix::splitAsBumpyMatrix(
#'     trans_info[, c("x", "y")], 
#'     row = trans_info$feature_name, col = trans_info$cell )
#' spe<- SpatialExperiment(
#'      assays = list(molecules = trans_mol),sample_id ="sample1" )
#' w_x <- c(min(floor(min(trans_info$x)),
#'          floor(min(clusters$x))),
#'       max(ceiling(max(trans_info$x)),
#'           ceiling(max(clusters$x))))
#' w_y <- c(min(floor(min(trans_info$y)),
#'           floor(min(clusters$y))),
#'       max(ceiling(max(trans_info$y)),
#'           ceiling(max(clusters$y))))
#' vecs_lst <- get_vectors(x=spe,sample_names=c("sample1"),
#'                     cluster_info = clusters,
#'                     bin_type = "square",
#'                     bin_param = c(20,20),
#'                     test_genes =c("gene_A1","gene_A2","gene_B1","gene_B2"),
#'                     w_x = w_x, w_y=w_y)
#' lasso_res <- lasso_markers(gene_mt=vecs_lst$gene_mt,
#'                         cluster_mt = vecs_lst$cluster_mt,
#'                         sample_names=c("sample1"),
#'                         keep_positive=TRUE,
#'                         background=NULL)
#' # the full result
#' full_result <- get_full_mg(lasso_res, coef_cutoff=0.05)
get_full_mg <- function(obj, coef_cutoff=0) {
    stopifnot(inherits(obj, "glm_mg_result"))
    t_res <-  obj$full_result
    t_res <- t_res[abs(t_res$glm_coef)>coef_cutoff,]
    t_res <- t_res[order(t_res$gene, t_res$p_value, na.last = TRUE), ]
    return(t_res)
}


