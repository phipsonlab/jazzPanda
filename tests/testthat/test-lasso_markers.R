
set.seed(12)
# simulate coordiantes for clusters
df_clA = data.frame(x = rnorm(n=100, mean=20, sd=5),
                    y = rnorm(n=100, mean=20, sd=5), cluster="A")
df_clB = data.frame(x = rnorm(n=100, mean=100, sd=5),
                    y = rnorm(n=100, mean=100, sd=5), cluster="B")

clusters = rbind(df_clA, df_clB)
clusters$sample="rep1"
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

w_x =  c(min(floor(min(trans_info$x)),
             floor(min(clusters$x))),
         max(ceiling(max(trans_info$x)),
             ceiling(max(clusters$x))))
w_y =  c(min(floor(min(trans_info$y)),
             floor(min(clusters$y))),
         max(ceiling(max(trans_info$y)),
             ceiling(max(clusters$y))))

data = list(trans_info = trans_info)
vecs_lst = get_vectors(data_lst=list(rep1=data), cluster_info = clusters,
                       bin_type = "square",
                       bin_param = c(20,20),
                       all_genes =c("gene_A1","gene_A2","gene_B1","gene_B2"),
                       w_x = w_x, w_y=w_y)
set.seed(100)
lasso_res1 = lasso_markers(gene_mt=vecs_lst$gene_mt,
                                cluster_mt = vecs_lst$cluster_mt,
                                sample_names=c("rep1"),
                                keep_positive=TRUE,
                                coef_cutoff=0.05,
                                background=NULL)
########################################################################
test_that("lasso- one sample case", {
  expect_equal(lasso_res1$lasso_top_result$top_cluster, c("A","A","B","B"))
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
                  coef_cutoff=0.05,
                  background=NULL))
})

########################################################################
set.seed(100)
test_that("Return error when some genes non-zero entry < n_fold", {
  expect_error(
    lasso_markers(gene_mt= vecs_lst$gene_mt,
                  cluster_mt =vecs_lst$cluster_mt,
                  sample_names=c("rep1"),
                  keep_positive=TRUE,
                  coef_cutoff=0.05,n_fold = 20,
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
                  coef_cutoff=0.05,
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
                  coef_cutoff=0.05,
                  background=NULL))
})

########################################################################
# test the randomness in cv for lasso
set.seed(989)
seed_lst = sample(1:9999, size = 100, replace = FALSE)
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
                               coef_cutoff=0.05,
                               background=NULL)
    all_res[(4*(curr_id-1)+1):(curr_id*4), ] = curr_res$lasso_top_result
    
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
                         coef_cutoff=0.05,
                         background=NULL)
one_gene_top = one_gene_lst$lasso_top_result
one_gene_full = one_gene_lst$lasso_full_result

test_that("lasso cross validation penalty - randomness", {
    expect_equal(unique(one_gene_top[one_gene_top$gene == "gene_A",
                                     "top_cluster"]), "A")
    expect_equal(unique(one_gene_full[one_gene_full$gene == "gene_A",
                                     "cluster"]), "A")
})

