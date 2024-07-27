#' helper function to check the inputs passed to marker detection function
#'
#'
#' @param gene_mt A matrix contains the transcript count in each grid.
#' Each row refers to a grid, and each column refers to a gene.
#' The column names must be specified and refer to the genes. This can be the
#' output from the function \code{\link{get_vectors}}.
#' @param cluster_mt A matrix contains the number of cells in a specific
#' cluster in each grid. Each row refers to a grid, and each column
#' refers to a cluster. The column names must be specified and refer to the
#' clusters. Please do not assign integers as column names.
#' This can be the output from the function \code{\link{get_vectors}}.
#' @param sample_names  A vector specifying the names for the samples.
#' @param n_fold Optional. A positive number giving the number of folds used
#' for cross validation. This parameter will pass to \code{\link{cv.glmnet}}
#' to calculate a penalty term for every gene.
#' @param background Optional. A matrix providing the
#' background information. Each row refers to a grid, and each column refers to
#' one category of background information. Number of rows must equal to the
#' number of rows in \code{gene_mt} and \code{cluster_mt}.
#' Can be obtained by only providing coordinates matrices \code{cluster_info}.
#' to function \code{get_vectors}.

#' @return a list of two matrices with the following components
#' \item{\code{n_clusters} }{ Number of clusters}
#' \item{\code{cluster_names} }{a vector of strings giving the name
#' of the clusters }

check_valid_input<- function(gene_mt,cluster_mt,sample_names, n_fold=10,
                            background=NULL){
    if (is.null(colnames(gene_mt)) ==TRUE){
        stop("Please provide valid column names for parameter gene_mt ") }
    if (is.null(colnames(cluster_mt)) ==TRUE){
        stop("Please provide valid column names for parameter cluster_mt ") }
    tryCatch({
        ot <- as.numeric(colnames(cluster_mt))
        if (all(is.na(ot))== FALSE){
            stop("Column names for cluster_mt contain integers. ") } })
    if (nrow(gene_mt) != nrow(cluster_mt)){
        stop("Number of rows of gene_mt and cluster_mt does not match.") }
    # every gene should have non-zero variance
    if ( any(colSums(gene_mt) == 0)){stop("Some genes have zero variance.")}
    # every cluster should have non-zero variance
    if (any(colSums(cluster_mt)==0)){stop("Some clusters have zero variance.")}
    # every gene should have at least n_fold non-zero entries
    if ( any(colSums(gene_mt !=0) < n_fold)){
        stop("Each gene should have at least n_fold non-zero entries") }
    # every cluster should have at least n_fold non-zero entries
    if ( any(colSums(cluster_mt !=0) < n_fold)){
        stop("Each cluster should have at least n_fold non-zero entries") }
    if (is.null(background) == FALSE){
        if (nrow(background) != nrow(gene_mt)){
            stop("Parameter background is in a wrong dimension. ") }
        if (length(colnames(background)) != ncol(background)){
            stop("Parameter background does not have valid column names. ") }
        if (length(intersect(colnames(background),colnames(cluster_mt))) != 0){
            stop("Parameter background and cluster_mt
                    have overlapped column names") } }
    if (length(sample_names) == 0 ){
        stop("Please provide a list of strings for parameter sample_names") }
    if (length(sample_names) == 1){
        n_clusters <- ncol(cluster_mt)
        cluster_names <- colnames(cluster_mt)
    }else{
        # check samples names are the same in cluster_mt for multiple samples
        if (length(setdiff(sample_names,colnames(cluster_mt))) != 0){
            stop("Inconsistent sample_names and sample names
                    from cluster_mt") }
        n_clusters <- ncol(cluster_mt) - length(sample_names)
        cluster_names <- setdiff(colnames(cluster_mt), sample_names) }
    return (list(n_clusters=n_clusters,cluster_names=cluster_names))
}

#' help function to get lasso coefficient for every cluster for a given model
#'
#'
#' @param i_gene Name of the current gene
#' @param gene_mt A matrix contains the transcript count in each grid.
#' Each row refers to a grid, and each column refers to a gene. The column
#' names must be specified and refer to the genes. This can be the
#' output from the function \code{\link{get_vectors}}.

#' @param vec_cluster  A matrix of the spatial vectors for clusters.

#' @param cluster_names  A vector of strings giving the name of clusters

#' @param n_fold Optional. A positive number giving the number of folds used
#' for cross validation. This parameter will pass to \code{\link{cv.glmnet}}
#' to calculate a penalty term for every gene.
#' @param n_samples A positive number giving the number samples
#' @param sample_names  A vector specifying the names for the sample
#'
#' @return a list of two matrices with the following components
#' \item{\code{coef_df} }{ A matrix giving the lasso coefficient of
#' each cluster}
#' \item{\code{lambda.1se} }{the lambda.1se value of best fitted model }
#'
get_lasso_coef <- function(i_gene, gene_mt,vec_cluster,cluster_names,n_fold=10,
                        n_samples, sample_names){

    data <- as.data.frame(matrix(gene_mt[,i_gene],ncol=1))
    colnames(data) <- "gene"
    # create stratified sampling
    ## add dummy noise to every gene vector value
    data$mod_count <- data$gene+runif(length(data$gene),
                                    min = 0.001, max = 0.01)
    ## get quantile of modified gene vector
    data$qt <- with(data, cut(mod_count,
                    breaks=quantile(mod_count,probs=seq(0,1, by=0.2)),
                    include.lowest=TRUE, labels=FALSE))

    folds <- createFolds(data$qt, k=n_fold,list = TRUE, returnTrain = FALSE)

    ## get the fold_id for every bin
    fold_id <- rep(NA, nrow(data))
    for(fd in seq_along(folds)){
        fold_id[folds[[fd]]] <- fd
    }
    if (n_samples == 1){
        lasso_reg <- cv.glmnet(x=vec_cluster[,cluster_names],
                        y=gene_mt[,i_gene],
                        family = "gaussian",
                        alpha =1, standardize = TRUE, intercept=FALSE,
                        keep=TRUE,
                        foldid = fold_id)

        best_model <- glmnet(x=vec_cluster[,cluster_names],
                            y=gene_mt[,i_gene],family = "gaussian",
                            alpha = 1, standardize = TRUE,
                            lambda = lasso_reg$lambda.min,intercept=FALSE)

    }else{
        lasso_reg <- cv.glmnet(x=vec_cluster[,c(cluster_names, sample_names)],
                            y=gene_mt[,i_gene],
                            family = "gaussian",
                            alpha =1, standardize = TRUE, intercept=FALSE,
                            keep=TRUE,
                            foldid = fold_id)

        best_model <- glmnet(x=vec_cluster[,c(cluster_names, sample_names)],
                            y=gene_mt[,i_gene],family = "gaussian",
                            alpha = 1, standardize = TRUE,
                            lambda = lasso_reg$lambda.min,intercept=FALSE)
    }
    coef_df <- as.matrix(coef(best_model))
    return (list(coef_df =coef_df,lambda.1se = lasso_reg$lambda.1se))
}

#' Find marker genes with spatial coordinates
#'
#' @description
#' This function will find the most spatially relevant cluster label
#' for each gene.
#'
#' @details
#' This function will take the converted gene and cluster vectors from function
#' \code{\link{get_vectors}}, and return the most relevant cluster label for
#' each gene. If there are multiple samples in the dataset, this function
#' will find shared markers across different samples by including additional
#' sample vectors in the input \code{cluster_mt}.
#'
#' This function treats all input cluster vectors as features, and create
#' a penalized linear model for one gene vector with lasso regularization.
#' Clusters with non-zero coefficient will be selected, and these clusters will
#' be used to formulate a generalised linear model for this gene vector.
#'
#' \itemize{
#' \item{If the input \code{keep_positive} is TRUE, the clusters
#' with positive coefficient and significant p-value will be saved in the
#' output matrix \code{lasso_full_result}. The cluster with a
#' positive coefficient and the minimum p-value will be regarded as the most
#' relevant cluster to this gene and be saved in the output matrix
#' \code{lasso_result}.}
#'
#' \item{If the input \code{keep_positive} is FALSE, the clusters with
#' negative coefficient and significant p-value will be saved in the
#' output matrix \code{lasso_full_result}. The cluster with a negative
#' coefficient and the minimum p-value will be regarded as the most relevant
#' cluster to this gene and be saved in the output matrix \code{lasso_result}.}
#' }
#'
#' If there is no clusters with significant p-value, the a string "NoSig" will
#' be returned for this gene.
#'
#' The parameter \code{background} can be used to capture unwanted noise
#' pattern in the dataset. For example, we can include negative control
#' genes as a background cluster in the model. If the most relevant cluster
#' selected by one gene matches the background "clusters",
#' we will return "NoSig" for this gene.
#'
#' @param gene_mt A matrix contains the transcript count in each grid.
#' Each row refers to a grid, and each column refers to a gene.
#' The column names must be specified and refer to the genes. This can be the
#' output from the function \code{\link{get_vectors}}.
#' @param cluster_mt A matrix contains the number of cells in a specific
#' cluster in each grid. Each row refers to a grid, and each column refers
#' to a cluster. The column names must be specified and refer to the clusters.
#' Please do not assign integers as column names.
#' This can be the output from the function \code{\link{get_vectors}}.
#' @param sample_names  A vector specifying the names for the samples.
#' @param keep_positive A logical flag indicating whether to return positively
#' correlated clusters or not.
#' @param coef_cutoff A positive number giving the coefficient cutoff value.
#' Genes whose top cluster showing a coefficient vlaue smaller than the cutoff
#' will be . Default is 0.05.
#' @param background Optional. A matrix providing the
#' background information. Each row refers to a grid, and each column refers to
#' one category of background information. Number of rows must equal to the
#' number of rows in \code{gene_mt} and \code{cluster_mt}.
#' Can be obtained by only providing coordinates matrices \code{cluster_info}.
#' to function \code{get_vectors}.
#' @param n_fold Optional. A positive number giving the number of folds used
#' for cross validation. This parameter will pass to \code{\link{cv.glmnet}}
#' to calculate a penalty term for every gene.

#' @return a list of two matrices with the following components
#' \item{\code{lasso_top_result}  }{A matrix with detailed information for
#' each gene and the most relevant cluster label.
#' \itemize{
#' \item{\code{gene} Gene name}
#' \item{\code{top_cluster} The name of the most revelant cluster
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
#'
#' }
#' \item{\code{lasso_full_result}  }{A matrix with detailed information for
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
#' }}

#' @importFrom caret createFolds
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet glmnet
#' @importFrom stats glm
#' @importFrom stats coef
#' @importFrom stats cor
#' @importFrom stats runif
#' @export
#'
#' @seealso \code{\link{get_vectors}}
#' @examples
#'
#' set.seed(100)
#' #  simulate coordinates for clusters
#' df_clA = data.frame(x = rnorm(n=100, mean=20, sd=5),
#'                  y = rnorm(n=100, mean=20, sd=5), cluster="A")
#' df_clB = data.frame(x = rnorm(n=100, mean=100, sd=5),
#'                 y = rnorm(n=100, mean=100, sd=5), cluster="B")
#'
#' clusters = rbind(df_clA, df_clB)
#' clusters$sample="rep1"
#'
#' # simulate coordinates for genes
#' trans_info = data.frame(rbind(cbind(x = rnorm(n=100, mean=20,sd=5),
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
#' trans_info$x=as.numeric(trans_info$x)
#' trans_info$y=as.numeric(trans_info$y)
#' w_x =  c(min(floor(min(trans_info$x)),
#'          floor(min(clusters$x))),
#'       max(ceiling(max(trans_info$x)),
#'           ceiling(max(clusters$x))))
#' w_y =  c(min(floor(min(trans_info$y)),
#'           floor(min(clusters$y))),
#'       max(ceiling(max(trans_info$y)),
#'           ceiling(max(clusters$y))))
#' data = list(trans_info = trans_info)
#' vecs_lst = get_vectors(data_lst=list(rep1=data), cluster_info = clusters,
#'                     bin_type = "square",
#'                     bin_param = c(20,20),
#'                     all_genes =c("gene_A1","gene_A2","gene_B1","gene_B2"),
#'                     w_x = w_x, w_y=w_y)
#' lasso_res = lasso_markers(gene_mt=vecs_lst$gene_mt,
#'                         cluster_mt = vecs_lst$cluster_mt,
#'                         sample_names=c("rep1"),
#'                         keep_positive=TRUE,
#'                         coef_cutoff=0.05,
#'                         background=NULL)
#'
#'
lasso_markers<- function(gene_mt,cluster_mt,sample_names,
                        keep_positive=TRUE, coef_cutoff=0.05,
                        background=NULL,n_fold=10){
    input_res <- check_valid_input(gene_mt=gene_mt, cluster_mt=cluster_mt,
                                sample_names=sample_names,n_fold=n_fold,
                                background=background)
    n_clusters <- input_res$n_clusters
    cluster_names <- input_res$cluster_names
    all_genes <- colnames(gene_mt)
    n_samples <- length(sample_names)
    n_genes <- length(all_genes)
    gene_mt <- as.matrix(gene_mt)
    nc_names <- NULL
    vec_cluster <- as.matrix(cluster_mt)
    # prepare background vectors
    if (is.null(background) == FALSE){
        nc_names<-colnames(background)
        vec_cluster<-as.matrix(cbind(vec_cluster, background))}
    res_df_all<-data.frame(matrix("", ncol=7))
    colnames(res_df_all)<-c("i_gene","cluster","Estimate","Pr(>|t|)","corr_p",
                                "max_gg_corr","max_gc_corr")
    lambda_vals <- matrix(0, nrow=n_genes, ncol=2)
    res_df<-as.data.frame(cbind(all_genes, matrix(0,nrow=n_genes,ncol=5)))
    colnames(res_df)<-c("gene","top1","glm_coef","pearson","max_gg_corr",
                            "max_gc_corr")
    row.names(res_df) <- all_genes
    res_df$gene <- as.character(res_df$gene)
    res_df$top1 <- as.character(res_df$top1)
    res_df$glm_coef <- as.numeric(res_df$glm_coef)
    res_df$max_gg_corr <- as.numeric(res_df$max_gg_corr)
    res_df$max_gc_corr <- as.numeric(res_df$max_gc_corr)
    res_df$pearson <- as.numeric(res_df$pearson)
    gene_cluster_m <- cor(x=gene_mt, y=vec_cluster, method = "pearson")
    gene_gene_corr <- cor(gene_mt,gene_mt,method = "pearson")
    if (ncol(gene_mt) > 1 ){ diag(gene_gene_corr)<- NA }
    gg_corr<-as.data.frame(cbind(apply(gene_gene_corr, MARGIN=1,
                        FUN = max, na.rm=TRUE),
                        apply(gene_cluster_m, MARGIN=1,FUN = max, na.rm=TRUE)))
    colnames(gg_corr) <- c("max_gg_corr","max_gc_corr")
    gg_corr$gene <- row.names(gg_corr)
    for (i in seq_len(ncol(gene_mt))){
        i_gene <- colnames(gene_mt)[i]
        gene_indx <- gg_corr$gene==i_gene
        res_df[i_gene,"max_gg_corr"]<-gg_corr[gene_indx,"max_gg_corr"]
        res_df[i_gene,"max_gc_corr"]<-gg_corr[gene_indx,"max_gc_corr"]
        lasso_coef_lst <- get_lasso_coef(i_gene=i_gene, gene_mt=gene_mt,
                vec_cluster=vec_cluster, cluster_names=cluster_names,
                n_fold=n_fold, n_samples=n_samples, sample_names=sample_names)
        lambda_vals[i,1] <- i_gene
        lambda_vals[i,2] <- lasso_coef_lst$lambda.1se
        coef_df <- lasso_coef_lst$coef_df
        if (length(coef_df[coef_df!=0]) ==0){
            res_df_all <- rbind(res_df_all,
                                cbind(i_gene,cluster="NoSig", Estimate=NA,
                                    `Pr(>|t|)`=NA,corr_p=NA,
                                    max_gg_corr=res_df[i_gene,"max_gg_corr"],
                                    max_gc_corr=res_df[i_gene,"max_gc_corr"]))
            res_df[i_gene,2] <- "NoSig"
            res_df[i_gene,3] <- 0
            res_df[i_gene,4] <- 0
        }else{
            coef_df <- coef_df[coef_df!=0,,drop=FALSE]
            if (n_samples>1){
                # multiple samples check number of significant clusterS
                if ( is.null(background) == TRUE){
                    # include sample variables as features
                    sig_features <- union(row.names(coef_df),sample_names)
                    glm_data <- data.frame("gene"=gene_mt[,i_gene],
                                        vec_cluster[,sig_features])
                }else{
                    # include sample variables and background info as features
                    sig_features <- union(row.names(coef_df),sample_names)
                    glm_data <- data.frame("gene"=gene_mt[,i_gene],
                                    vec_cluster[,c(sig_features, nc_names)])}
                coef_mod<-as.matrix(summary(glm(gene~.+0,
                            data=glm_data,family = "gaussian"))$coefficients)
                # positive coefficient
                if (keep_positive == TRUE){
                    target_clusters <- coef_mod[coef_mod[,"Pr(>|t|)"]<= 0.05 &
                                coef_mod[,"Estimate"]>0,c(1,4),drop=FALSE]
                }else{
                    target_clusters<-coef_mod[coef_mod[,"Pr(>|t|)"]<= 0.05 &
                                coef_mod[,"Estimate"]<0, c(1,4),drop=FALSE] }
                ccluster_count<-length(setdiff(row.names(target_clusters),
                                                sample_names))
        }else{
                # one-sample, check the number of significant clusters
                if (is.null(background) == TRUE){
                    glm_data<-data.frame("gene"=gene_mt[,i_gene],
                                vec_cluster[,row.names(coef_df),drop=FALSE])
                }else{
                    # include background info as features
                    glm_data<-data.frame("gene"=gene_mt[,i_gene],
                            vec_cluster[,c(row.names(coef_df), nc_names),
                                            drop=FALSE]) }
                coef_mod<-as.matrix(summary(glm(gene~.+0, data=glm_data,
                                            family = "gaussian"))$coefficients)
                if (keep_positive == TRUE){
                    # positive coefficient
                    target_clusters<-coef_mod[coef_mod[,"Pr(>|t|)"]<= 0.05 &
                                    coef_mod[,"Estimate"]>0, c(1,4),drop=FALSE]
                }else{
                    # negative coefficient
                    target_clusters<-coef_mod[coef_mod[,"Pr(>|t|)"]<= 0.05 &
                                    coef_mod[,"Estimate"]<0,c(1,4),drop=FALSE]}
                ccluster_count<-length(row.names(target_clusters)) }
            # set gene as NA if no clusters are significant
            if ((nrow(target_clusters) == 0) | (ccluster_count==0) ) {
                res_df_all<-rbind(res_df_all,
                                cbind(i_gene, cluster="NoSig", Estimate=NA,
                                    `Pr(>|t|)`=NA, corr_p=NA,
                                    max_gg_corr=res_df[i_gene,"max_gg_corr"],
                                    max_gc_corr=res_df[i_gene,"max_gc_corr"]))
                res_df[i_gene,2] <- "NoSig"
                res_df[i_gene,3] <- 0
                res_df[i_gene,4] <- 0
            }else{
                sig_clusters<-as.matrix(target_clusters)
                sig_clusters<-sig_clusters[order(sig_clusters[,2],
                                    -abs(sig_clusters[,1])),,drop=FALSE]
                corr_p<-as.vector(cor(gene_mt[,i_gene],
                    vec_cluster[,row.names(sig_clusters)],method = "pearson"))
                res_df_all<-rbind(res_df_all, cbind(i_gene,
                        cluster=row.names(sig_clusters), sig_clusters,corr_p,
                        max_gg_corr=res_df[i_gene,"max_gg_corr"],
                        max_gc_corr=res_df[i_gene,"max_gc_corr"]))
                sigc_nms <- row.names(sig_clusters)
                res_df[i_gene,2]<-setdiff(sigc_nms,sample_names)[1]
                res_df[i_gene,3]<-sig_clusters[res_df[i_gene,2],"Estimate"]
                res_df[i_gene,4]<-cor(gene_mt[,i_gene],
                        vec_cluster[,res_df[i_gene,2]],method = "pearson")} } }
    res_df$top_cluster <- res_df$top1
    res_df[abs(res_df$glm_coef)<= coef_cutoff,"top_cluster"] <- "NoSig"
    res_df[abs(res_df$glm_coef)<= coef_cutoff,"glm_coef"] <- 0
    res_df[abs(res_df$glm_coef)<= coef_cutoff,"pearson"] <- 0
    res_df<-res_df[, c("gene","top_cluster","glm_coef","pearson",
                        "max_gg_corr","max_gc_corr")]
    colnames(res_df_all)<-c("gene","cluster","glm_coef","p_value","pearson",
                        "max_gg_corr","max_gc_corr")
    if (is.null(background) == FALSE){
        nc_rows_id <- which(res_df$top_cluster %in% nc_names)
        res_df[nc_rows_id, "top_cluster"] <- "NoSig"
        res_df[nc_rows_id, "glm_coef"] <- 0
        res_df[nc_rows_id, "pearson"] <- 0 }
    res_df_all <- res_df_all[2:nrow(res_df_all),]
    #res_df_all <- res_df_all[res_df_all$cluster %in% cluster_names,]
    row.names(res_df_all) <- NULL
    res_df_all$glm_coef  <- as.numeric(res_df_all$glm_coef)
    res_df_all$p_value  <- as.numeric(res_df_all$p_value)
    res_df_all$pearson  <- as.numeric(res_df_all$pearson)
    res_df_all$glm_coef  <- as.numeric(res_df_all$glm_coef)
    to_num_cols<-c("glm_coef","p_value","pearson","max_gg_corr","max_gc_corr")
    res_df_all[, to_num_cols] <- lapply(res_df_all[to_num_cols], as.numeric)
    res_df_all <- res_df_all[order(res_df_all$gene, res_df_all$p_value,
                                na.last = TRUE), ]
    return (list(lasso_top_result = data.frame(res_df),
                lasso_full_result = data.frame(res_df_all)))
}
