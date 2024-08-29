set.seed(100)
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

rep1 = list(trans_info = trans_info)


test_that("Invaid input",{
    expect_error(compute_permp(data=rep1,
                               cluster_info=clusters,
                               perm.size=100,
                               bin_type="square",
                               bin_param=c(2),
                               all_genes=unique(trans_info$feature_name),
                               correlation_method = "pearson",
                               n.cores=2,
                               correction_method="BH",
                               w_x=w_x ,
                               w_y=w_y))
    expect_error(compute_permp(data=rep1,
                               cluster_info=clusters,
                               perm.size=100,
                               bin_type="hexagon",
                               bin_param=c(2,2),
                               all_genes=unique(trans_info$feature_name),
                               correlation_method = "pearson",
                               n.cores=2,
                               correction_method="BH",
                               w_x=w_x ,
                               w_y=w_y))
    expect_error(compute_permp(data=rep1,
                               cluster_info=clusters,
                               perm.size=100,
                               bin_type="circle",
                               bin_param=c(2,2),
                               all_genes=unique(trans_info$feature_name),
                               correlation_method = "pearson",
                               n.cores=2,
                               correction_method="BH",
                               w_x=w_x ,
                               w_y=w_y))
    
})
set.seed(100)
perm_p_lst = compute_permp(data=rep1,
                       cluster_info=clusters,
                       perm.size=100,
                       bin_type="square",
                       bin_param=c(2,2),
                       all_genes=unique(trans_info$feature_name),
                       correlation_method = "pearson",
                       n.cores=2,
                       correction_method="BH",
                       w_x=w_x ,
                       w_y=w_y)
perm_p_s = compute_permp(data=rep1,
                           cluster_info=clusters,
                           perm.size=100,
                           bin_type="square",
                           bin_param=c(2,2),
                           all_genes=unique(trans_info$feature_name),
                           correlation_method = "pearson",
                           n.cores=1,
                           correction_method="BH",
                           w_x=w_x ,
                           w_y=w_y)
test_that("Test permutation result - output dimension matches", {
  expect_equal(length(perm_p_lst), 4)
  expect_equal(dim(perm_p_lst$perm.arrays), c(4,2,100))
  expect_equal(dim(perm_p_lst$obs.stat), c(4,2))
  expect_equal(dim(perm_p_lst$perm.pval.adj), c(4,2))
  expect_equal(dim(perm_p_lst$perm.pval), c(4,2))
  expect_equal(names(perm_p_lst),
               c("obs.stat", "perm.arrays", "perm.pval", "perm.pval.adj"))
  

})

test_that("Test permutation result - observed stat matches", {
  expect_equal(as.vector(perm_p_lst$obs.stat),
               c(1,1, -1/3, -1/3,-1/3,-1/3,1,1))
})

perm_p_s = compute_permp(data=rep1,
                         cluster_info=clusters,
                         perm.size=100,
                         bin_type="square",
                         bin_param=c(2,2),
                         all_genes=unique(trans_info$feature_name),
                         correlation_method = "pearson",
                         n.cores=1,
                         correction_method="BH",
                         w_x=w_x ,
                         w_y=w_y)
test_that("Test permutation result - sequential calculation works", {
    expect_equal(length(perm_p_s), 4)
    expect_equal(dim(perm_p_s$perm.arrays), c(4,2,100))
    expect_equal(dim(perm_p_s$obs.stat), c(4,2))
    expect_equal(dim(perm_p_s$perm.pval.adj), c(4,2))
    expect_equal(dim(perm_p_s$perm.pval), c(4,2))
    expect_equal(names(perm_p_s),
                 c("obs.stat", "perm.arrays", "perm.pval", "perm.pval.adj"))
    expect_equal(as.vector(perm_p_s$obs.stat),
                 as.vector(perm_p_lst$obs.stat))
    
    
})

#############################################################################
set.seed(100)
perm_hex_lst = compute_permp(data=rep1,
                           cluster_info=clusters,
                           perm.size=100,
                           bin_type="hexagon",
                           bin_param=c(5),
                           all_genes=unique(trans_info$feature_name),
                           correlation_method = "pearson",
                           n.cores=2,
                           correction_method="BH",
                           w_x=w_x ,
                           w_y=w_y)
test_that("Test permutation result - output dimension matches", {
    expect_equal(length(perm_hex_lst), 4)
    expect_equal(dim(perm_hex_lst$perm.arrays), c(4,2,100))
    expect_equal(dim(perm_hex_lst$obs.stat), c(4,2))
    expect_equal(dim(perm_hex_lst$perm.pval.adj), c(4,2))
    expect_equal(dim(perm_hex_lst$perm.pval), c(4,2))
    expect_equal(names(perm_hex_lst),
                 c("obs.stat", "perm.arrays", "perm.pval", "perm.pval.adj"))
    
})

