library(SpatialExperiment)
set.seed(12)
# simulate coordiantes for clusters
df_clA = data.frame(x = rnorm(n=100, mean=20, sd=5),
                    y = rnorm(n=100, mean=20, sd=5), cluster="A")
df_clB = data.frame(x = rnorm(n=100, mean=100, sd=5),
                    y = rnorm(n=100, mean=100, sd=5), cluster="B")

clusters = rbind(df_clA, df_clB)
clusters$sample="sample1"
# simulate coordiantes for genes
trans_info = data.frame(rbind(cbind(x = rnorm(n=100, mean=20, sd=5),
                                    y = rnorm(n=100, mean=20, sd=5),
                                    feature_name="gene_A1"),
                              cbind(x = rnorm(n=100, mean=20, sd=5),
                                    y = rnorm(n=100, mean=20, sd=5),
                                    feature_name="gene_A2"),
                              cbind(x = rnorm(n=100, mean=100, sd=5),
                                    y = rnorm(n=100, mean=100, sd=5),
                                    feature_name="gene_B1"),
                              cbind(x = rnorm(n=100, mean=100, sd=5),
                                    y = rnorm(n=100, mean=100, sd=5),
                                    feature_name="gene_B2")))
trans_info$x=as.numeric(trans_info$x)
trans_info$y=as.numeric(trans_info$y)
trans_info$cell_id = paste("cell", 1:nrow(trans_info), sep="")
# create SpatialExperiment object
trans_mol <- BumpyMatrix::splitAsBumpyMatrix(
    trans_info[, c("x", "y")], 
    row = trans_info$feature_name, col = trans_info$cell_id )

spe<- SpatialExperiment(
    assays = list(molecules = trans_mol),sample_id ="sample1" )


vecs_lst = get_vectors(x=spe, cluster_info = clusters,
                       bin_type = "square",sample_names = "sample1",
                       bin_param = c(20,20),
                       test_genes =c("gene_A1","gene_A2","gene_B1","gene_B2"))
#### background

background_sv = create_genesets(x =spe,sample_names = "sample1",
                                name_lst=list(dummy_W=c("gene_A1","gene_B1")),
                                bin_type="square",
                                bin_param = c(20,20),cluster_info = NULL)
set.seed(100)
lasso_res1 = lasso_markers(gene_mt=vecs_lst$gene_mt,
                           cluster_mt = vecs_lst$cluster_mt,
                           sample_names=c("rep1"),
                           keep_positive=TRUE,
                           background=NULL)
set.seed(100)
lasso_res0 = lasso_markers(gene_mt=vecs_lst$gene_mt,
                           cluster_mt = vecs_lst$cluster_mt,
                           sample_names=c("rep1"),
                           keep_positive=TRUE,
                           background=NULL)

set.seed(100)
lasso_res_neg = lasso_markers(gene_mt=vecs_lst$gene_mt*(-1),
                              cluster_mt = vecs_lst$cluster_mt,
                              sample_names=c("rep1"),
                              keep_positive=FALSE,
                              background=NULL)
set.seed(100)
lasso_res_background = lasso_markers(gene_mt=vecs_lst$gene_mt*(-1),
                              cluster_mt = vecs_lst$cluster_mt,
                              sample_names=c("rep1"),
                              keep_positive=FALSE,
                              background=background_sv)
invlid_vec = vecs_lst
colnames(invlid_vec$cluster_mt)=NULL
colnames(invlid_vec$gene_mt)=NULL
# clusters are named as integers 
invlid_vec_int = vecs_lst
colnames(invlid_vec_int$cluster_mt)=1:ncol(invlid_vec_int$cluster_mt)
# background has NULL names
invalid_back =background_sv
colnames(invalid_back)=NULL
# background matrix has overlaps with cluster_mt
invalid_back_overlapped =vecs_lst$cluster_mt
# cluster_mt with < nfold nonzeros
invalid_mt_le_fold = cbind(vecs_lst$cluster_mt,
                           invalid_col=rep(0, nrow(vecs_lst$cluster_mt)))
invalid_mt_le_fold[1:8,ncol(invalid_mt_le_fold)] = 1:8
test_that("Invalid input", {
    expect_error(lasso_markers(gene_mt=vecs_lst$gene_mt[1:2, ],
                               cluster_mt = vecs_lst$cluster_mt,
                               sample_names=NULL,
                               keep_positive=FALSE,
                               background=NULL))
    expect_error(lasso_markers(gene_mt=vecs_lst$gene_mt,
                               cluster_mt = vecs_lst$cluster_mt,
                               sample_names=NULL,
                               keep_positive=FALSE,
                               background=NULL))
    expect_error(lasso_markers(gene_mt=vecs_lst$gene_mt,
                               cluster_mt = invlid_vec$cluster_mt,
                               sample_names=c("rep1"),
                               keep_positive=FALSE,
                               background=background_sv))
    expect_error(lasso_markers(gene_mt=invlid_vec$gene_mt,
                               cluster_mt = vecs_lst$cluster_mt,
                               sample_names=c("rep1"),
                               keep_positive=FALSE,
                               background=background_sv))
    expect_error(lasso_markers(gene_mt=vecs_lst$gene_mt,
                               cluster_mt = vecs_lst$cluster_mt,
                               sample_names=c("rep1"),
                               keep_positive=FALSE,
                               background=invalid_back))
    expect_error(lasso_markers(gene_mt=vecs_lst$gene_mt,
                               cluster_mt = vecs_lst$cluster_mt,
                               sample_names=c("rep1"),
                               keep_positive=FALSE,
                               background=invalid_back_overlapped))
    expect_error(lasso_markers(gene_mt=vecs_lst$gene_mt,
                               cluster_mt = vecs_lst$cluster_mt,
                               sample_names=c("rep1"),
                               keep_positive=FALSE,
                               background=background_sv[1:10,]))
    expect_error(lasso_markers(gene_mt=vecs_lst$gene_mt,
                               cluster_mt = invalid_mt_le_fold,
                               sample_names=c("rep1"),
                               keep_positive=FALSE,
                               background=NULL))
    expect_error(lasso_markers(gene_mt=vecs_lst$gene_mt,
                               cluster_mt = invlid_vec_int$cluster_mt,
                               sample_names=NULL,
                               keep_positive=FALSE,
                               background=NULL))
})

########################################################################
test_that("lasso- one sample case - dimension correct", {
    expect_equal(length(lasso_res1), 2)
})
top_res <- get_top_mg(lasso_res1,coef_cutoff = 0.05)
test_that("lasso- one sample case", {
    expect_equal(top_res$top_cluster, c("A","A","B","B"))
    
})
########################################################################
top_res <- get_top_mg(lasso_res_neg,coef_cutoff = 0.05)
full_res <- get_full_mg(lasso_res_neg,coef_cutoff = 0.05)
test_that("lasso- one sample case- 0 cutoff - dimension correct", {
    expect_equal(length(lasso_res_neg), 2)
})

test_that("lasso- one sample case - 0 cutoff", {
    expect_equal(top_res$top_cluster, c("A","A","B","B"))
    expect_equal(full_res$cluster, c("A","A","B","B"))
})
########################################################################
top_res <- get_top_mg(lasso_res0,coef_cutoff = 0.05)
full_res <- get_full_mg(lasso_res0,coef_cutoff = 0.05)

test_that("lasso- one sample case -negative correlation - dimension correct", {
    expect_equal(length(lasso_res0), 2)
})
test_that("lasso- one sample case -negative correlation", {
    expect_equal(top_res$top_cluster, c("A","A","B","B"))
    expect_equal(full_res$cluster, c("A","A","B","B"))
})
########################################################################
top_res <- get_top_mg(lasso_res_background,coef_cutoff = 0.05)
full_res <- get_full_mg(lasso_res_background,coef_cutoff = 0.05)

test_that("lasso- wih background", {
    expect_equal(length(lasso_res_background), 2)
    expect_equal(top_res$top_cluster, c("A","A","B","B"))
    expect_equal(unique(full_res$cluster), c("A","dummy_W","B"))
})
########################################################################
dummy_cl = vecs_lst$cluster_mt
colnames(dummy_cl) = c("1","B")
set.seed(100)
test_that("Return error when column names contain integers", {
    expect_error(
        lasso_markers(gene_mt=vecs_lst$gene_mt,
                      cluster_mt = dummy_cl,
                      sample_names=c("rep1"),
                      keep_positive=TRUE,
                      background=NULL))
})

########################################################################
set.seed(100)
test_that("Return error when some genes non-zero entry < n_fold", {
    expect_error(
        lasso_markers(gene_mt= vecs_lst$gene_mt,
                      cluster_mt =vecs_lst$cluster_mt,
                      sample_names=c("rep1"),
                      keep_positive=TRUE,n_fold = 20,
                      background=NULL))
})
########################################################################
set.seed(100)
dummy_cl = vecs_lst$cluster_mt
dummy_cl[,1] = 0
test_that("Return error when some clusters have zero variance", {
    expect_error(
        lasso_markers(gene_mt=vecs_lst$gene_mt,
                      cluster_mt = dummy_cl,
                      sample_names=c("rep1"),
                      keep_positive=TRUE,
                      background=NULL))
})

dummy_ge = vecs_lst$gene_mt
dummy_ge[,1] = 0
test_that("Return error when some genes have zero variance", {
    expect_error(
        lasso_markers(gene_mt=dummy_ge,
                      cluster_mt =vecs_lst$cluster_mt,
                      sample_names=c("rep1"),
                      keep_positive=TRUE,
                      background=NULL))
})

########################################################################
# test the randomness in cv for lasso
set.seed(989)
seed_lst = sample(1:9999, size = 5, replace = FALSE)
all_res = as.data.frame(matrix(0, nrow = 4*length(seed_lst), ncol = 6))
colnames(all_res)= c("gene","top_cluster","glm_coef","pearson",
                     "max_gg_corr","max_gc_corr")

for (curr_id in 1:length(seed_lst)){
    sed = seed_lst[curr_id]
    set.seed(sed)
    curr_res = lasso_markers(gene_mt=vecs_lst$gene_mt,
                             cluster_mt = vecs_lst$cluster_mt,
                             sample_names=c("rep1"),
                             keep_positive=TRUE,
                             background=NULL)
    top_res = get_top_mg(curr_res, coef_cutoff = 0.05)
    all_res[(4*(curr_id-1)+1):(curr_id*4), ] = top_res
    
}

test_that("lasso cross validation penalty - randomness", {
    expect_equal(unique(all_res[all_res$gene == "gene_A1","top_cluster"]), "A")
    expect_equal(unique(all_res[all_res$gene == "gene_A2","top_cluster"]), "A")
    expect_equal(unique(all_res[all_res$gene == "gene_B1","top_cluster"]), "B")
    expect_equal(unique(all_res[all_res$gene == "gene_B2","top_cluster"]), "B")
})


########################################################################
# test case with only one gene 
set.seed(989)
one_gene_mt  = vecs_lst$gene_mt[,1, drop=FALSE]
colnames(one_gene_mt) = "gene_A"
one_gene_lst = lasso_markers(gene_mt=one_gene_mt,
                             cluster_mt = vecs_lst$cluster_mt,
                             sample_names=c("rep1"),
                             keep_positive=TRUE,
                             background=NULL)
one_gene_top = get_top_mg(one_gene_lst, coef_cutoff = 0.05)
one_gene_full = get_full_mg(one_gene_lst, coef_cutoff = 0.05)

test_that("lasso cross validation penalty - randomness", {
    expect_equal(unique(one_gene_top[one_gene_top$gene == "gene_A",
                                     "top_cluster"]), "A")
    expect_equal(unique(one_gene_full[one_gene_full$gene == "gene_A",
                                      "cluster"]), "A")
})

one_gene_res1 = lasso_markers(gene_mt=one_gene_mt,
                              cluster_mt = vecs_lst$cluster_mt,
                              sample_names=c("rep1"),
                              keep_positive=TRUE,
                              background=NULL)
one_gene_top = get_top_mg(one_gene_res1, coef_cutoff = 2)
one_gene_full =get_full_mg(one_gene_res1, coef_cutoff = 0)
test_that("high resolution results in 0 significant genes", {
    expect_equal(unique(one_gene_top[one_gene_top$gene == "gene_A",
                                     "top_cluster"]), "NoSig")
    expect_equal(unique(one_gene_full[one_gene_full$gene == "gene_A",
                                      "cluster"]), "A")
})


###############################################################################    
# test two sample scenario
set.seed(12)
# simulate coordiantes for clusters
df_clA = data.frame(x = rnorm(n=200, mean=20, sd=5),
                    y = rnorm(n=200, mean=20, sd=5), cluster="A")
df_clB = data.frame(x = rnorm(n=200, mean=100, sd=5),
                    y = rnorm(n=200, mean=100, sd=5), cluster="B")

clusters = rbind(df_clA, df_clB)
clusters$sample=c(rep("sample1",times=100),rep("sample2",times=100))
# simulate coordiantes for genes
trans_info_sp2 = data.frame(rbind(cbind(x = rnorm(n=100, mean=20, sd=5),
                                    y = rnorm(n=100, mean=20, sd=5),
                                    feature_name="gene_A1"),
                              cbind(x = rnorm(n=100, mean=20, sd=5),
                                    y = rnorm(n=100, mean=20, sd=5),
                                    feature_name="gene_A2"),
                              cbind(x = rnorm(n=100, mean=100, sd=5),
                                    y = rnorm(n=100, mean=100, sd=5),
                                    feature_name="gene_B1"),
                              cbind(x = rnorm(n=100, mean=100, sd=5),
                                    y = rnorm(n=100, mean=100, sd=5),
                                    feature_name="gene_B2")))
trans_info_sp2$x=as.numeric(trans_info_sp2$x)
trans_info_sp2$y=as.numeric(trans_info_sp2$y)

trans_info_sp2$cell_id = paste("cell", 1:nrow(trans_info_sp2), sep="")
# create SpatialExperiment object
sp2_trans_mol <- BumpyMatrix::splitAsBumpyMatrix(
    trans_info_sp2[, c("x", "y")], 
    row = trans_info_sp2$feature_name, col = trans_info_sp2$cell_id )

spe_sp2<- SpatialExperiment(
    assays = list(molecules = sp2_trans_mol),sample_id ="sample2" )
# 
# w_x =  c(min(floor(min(trans_info$x)),floor(min(trans_info_sp2$x)),
#              floor(min(clusters$x))),
#          max(ceiling(max(trans_info$x)),ceiling(max(trans_info_sp2$x)),
#              ceiling(max(clusters$x))))
# w_y =  c(min(floor(min(trans_info$y)),floor(min(trans_info_sp2$y)),
#              floor(min(clusters$y))),
#          max(ceiling(max(trans_info$y)),ceiling(max(trans_info_sp2$y)),
#              ceiling(max(clusters$y))))

twosample = SingleCellExperiment::cbind(spe,spe_sp2)
vecs_lst = get_vectors(x=twosample, sample_names = c("sample1","sample2"),
                       cluster_info = clusters,
                       bin_type = "square",
                       bin_param = c(20,20),
                       test_genes =c("gene_A1","gene_A2","gene_B1","gene_B2"))
#### background

background_sv = create_genesets(x=twosample, sample_names = c("sample1","sample2"),
                                name_lst=list(dummy_W=c("gene_A1","gene_B1")),
                                bin_type="square",
                                bin_param = c(20,20), cluster_info = NULL)
set.seed(100)
lasso_res_background = lasso_markers(gene_mt=vecs_lst$gene_mt,
                                     cluster_mt = vecs_lst$cluster_mt,
                                     sample_names=c("sample1","sample2"),
                                     keep_positive=TRUE,
                                     background=background_sv)
top_res <- get_top_mg(lasso_res_background,coef_cutoff = 0)
full_res <- get_full_mg(lasso_res_background,coef_cutoff = 0)

test_that("lasso-two sample with background", {
    expect_equal(length(lasso_res_background), 2)
    expect_equal(top_res$top_cluster, c("A","A","B","B"))
    expect_equal(unique(full_res$cluster), c("A","dummy_W","B"))
    expect_error( lasso_markers(gene_mt=vecs_lst$gene_mt,
                                cluster_mt = vecs_lst$cluster_mt,
                                sample_names=c("sample3","sample4"),
                                keep_positive=TRUE,
                                background=NULL))
})


#######
background_sv_all = create_genesets(x=twosample, 
                                    sample_names = c("sample1","sample2"),
                                name_lst=list(dummy_A1=c("gene_A1"),
                                              dummy_B1=c("gene_B1"),
                                              dummy_A2=c("gene_A2"),
                                              dummy_B2=c("gene_B2")),
                                bin_type="square",
                                bin_param = c(20,20),
                                cluster_info = NULL)
lasso_nosig_bg = lasso_markers(gene_mt=vecs_lst$gene_mt,
                                     cluster_mt = vecs_lst$cluster_mt,
                                     sample_names=c("sample1","sample2"),
                                     keep_positive=TRUE,
                                     background=background_sv_all)
top_res_bg <- get_top_mg(lasso_nosig_bg,coef_cutoff = 0)
full_res_bg <- get_full_mg(lasso_nosig_bg,coef_cutoff = 0)

set.seed(989)
rand_cl_mt = matrix(runif(n=1600,min=0, max = 1),
                    ncol=2, 
                    nrow=nrow(vecs_lst$cluster_mt))
colnames(rand_cl_mt)=colnames(vecs_lst$cluster_mt)[1:2]
rand_cl_mt = as.data.frame(cbind(rand_cl_mt, vecs_lst$cluster_mt[, c(3:4)]))
lasso_nosig_cl = lasso_markers(gene_mt=vecs_lst$gene_mt,
                            cluster_mt = rand_cl_mt,
                            sample_names=c("sample1","sample2"),
                            keep_positive=TRUE,
                            background=NULL)
top_res_cl <- get_top_mg(lasso_nosig_cl,coef_cutoff = 0)
full_res_cl <- get_full_mg(lasso_nosig_cl,coef_cutoff = 0)

test_that("set gene as NA if no clusters are significant", {
    expect_equal(length(lasso_nosig_bg), 2)
    expect_equal(unique(top_res_bg$top_cluster), c("NoSig"))
    expect_equal(unique(full_res_bg$cluster), 
                 c("dummy_A1","A","dummy_A2","dummy_B1","B","dummy_B2"))
    
    expect_equal(length(lasso_nosig_cl), 2)
    expect_equal(unique(top_res_cl$top_cluster), c("NoSig"))
})

