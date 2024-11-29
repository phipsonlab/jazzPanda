library(SpatialExperiment)
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

trans_info$cell =  rep(paste("cell",1:200, sep=""), times=2)
mol <- BumpyMatrix::splitAsBumpyMatrix(
    trans_info[, c("x", "y")], 
    row = trans_info$feature_name, col = trans_info$cell )

spe_rep1 <- SpatialExperiment(
    assays = list(molecules = mol),sample_id ="rep1" )
spe_rep2 <- SpatialExperiment(
    assays = list(molecules = mol),sample_id ="rep2" )

invalid_spe_rep1 <-  SingleCellExperiment::cbind(spe_rep1,spe_rep2)

w_x =  c(min(floor(min(trans_info$x)),
             floor(min(clusters$x))),
         max(ceiling(max(trans_info$x)),
             ceiling(max(clusters$x))))
w_y =  c(min(floor(min(trans_info$y)),
             floor(min(clusters$y))),
         max(ceiling(max(trans_info$y)),
             ceiling(max(clusters$y))))

test_that("Invaid input",{
    expect_error(compute_permp(x=trans_info,
                               cluster_info=clusters,
                               perm.size=10,
                               bin_type="square",
                               bin_param=c(2),
                               test_genes=unique(trans_info$feature_name),
                               correlation_method = "pearson",
                               n_cores=2,
                               correction_method="BH",
                               w_x=w_x ,
                               w_y=w_y))
    
    expect_error(compute_permp(x=spe_rep1,
                               cluster_info=clusters,
                               perm.size=100,
                               bin_type="square",
                               bin_param=c(2),
                               test_genes=unique(trans_info$feature_name),
                               correlation_method = "pearson",
                               n_cores=2,
                               correction_method="BH",
                               w_x=w_x ,
                               w_y=w_y))
    expect_error(compute_permp(x=spe_rep1,
                               cluster_info=clusters,
                               perm.size=100,
                               bin_type="hexagon",
                               bin_param=c(2,2),
                               test_genes=unique(trans_info$feature_name),
                               correlation_method = "pearson",
                               n_cores=2,
                               correction_method="BH",
                               w_x=w_x ,
                               w_y=w_y))
    expect_error(compute_permp(x=spe_rep1,
                               cluster_info=clusters,
                               perm.size=100,
                               bin_type="circle",
                               bin_param=c(2,2),
                               test_genes=unique(trans_info$feature_name),
                               correlation_method = "pearson",
                               n_cores=2,
                               correction_method="BH",
                               w_x=w_x ,
                               w_y=w_y))
    
    expect_error(compute_permp(x=invalid_spe_rep1,
                               cluster_info=clusters,
                               perm.size=100,
                               bin_type="circle",
                               bin_param=c(2,2),
                               test_genes=unique(trans_info$feature_name),
                               correlation_method = "pearson",
                               n_cores=2,
                               correction_method="BH",
                               w_x=w_x ,
                               w_y=w_y))
})
set.seed(100)
perm_p_lst = compute_permp(x=spe_rep1,
                       cluster_info=clusters,
                       perm.size=10,
                       bin_type="square",
                       bin_param=c(2,2),
                       test_genes=unique(trans_info$feature_name),
                       correlation_method = "pearson",
                       n_cores=2,
                       correction_method="BH",
                       w_x=w_x ,
                       w_y=w_y)
perm_p_s = compute_permp(x=spe_rep1,
                           cluster_info=clusters,
                         perm.size=10,
                           bin_type="square",
                           bin_param=c(2,2),
                           test_genes=unique(trans_info$feature_name),
                           correlation_method = "pearson",
                           n_cores=1,
                           correction_method="BH",
                           w_x=w_x ,
                           w_y=w_y)
test_that("Test permutation result - output dimension matches", {
  expect_equal(length(perm_p_lst), 4)
  expect_equal(dim(perm_p_lst$perm.arrays), c(4,2,10))
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

perm_p_s = compute_permp(x=spe_rep1,
                         cluster_info=clusters,
                         perm.size=10,
                         bin_type="square",
                         bin_param=c(2,2),
                         test_genes=unique(trans_info$feature_name),
                         correlation_method = "pearson",
                         n_cores=1,
                         correction_method="BH",
                         w_x=w_x ,
                         w_y=w_y)
test_that("Test permutation result - sequential calculation works", {
    expect_equal(length(perm_p_s), 4)
    expect_equal(dim(perm_p_s$perm.arrays), c(4,2,10))
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
perm_hex_lst = compute_permp(x=spe_rep1,
                           cluster_info=clusters,
                           perm.size=10,
                           bin_type="hexagon",
                           bin_param=c(5),
                           test_genes=unique(trans_info$feature_name),
                           correlation_method = "pearson",
                           n_cores=2,
                           correction_method="BH",
                           w_x=w_x ,
                           w_y=w_y)
test_that("Test permutation result - output dimension matches", {
    expect_equal(length(perm_hex_lst), 4)
    expect_equal(dim(perm_hex_lst$perm.arrays), c(4,2,10))
    expect_equal(dim(perm_hex_lst$obs.stat), c(4,2))
    expect_equal(dim(perm_hex_lst$perm.pval.adj), c(4,2))
    expect_equal(dim(perm_hex_lst$perm.pval), c(4,2))
    expect_equal(names(perm_hex_lst),
                 c("obs.stat", "perm.arrays", "perm.pval", "perm.pval.adj"))
    
})

#############################################################################
# 
# generate gene vector from count matrix
cm <- data.frame(rbind("gene_A"=c(0,0,2,0,0,0,2),
                       "gene_B"=c(5,3,3,13,0,1,14),
                       "gene_C"=c(5,0,1,5,1,0,7),
                       "gene_D"=c(0,1,1,2,0,0,2)))
colnames(cm)= paste("cell_", 1:7, sep="")

# simulate coordiantes for clusters
clusters = data.frame(x = c(1, 2,20,21,22,23,24),
                      y = c(23, 24, 1,2,3,4,5), cluster="A")
clusters$sample="rep1"
clusters$cell_id= colnames(cm)
# simulate coordiantes for genes
w_x=c(0,25)
w_y=c(0,25)
# cell_1 = (1,0,0,0)
# cell_2 = (1,0,0,0)
# cell_3 = (0,0,0,1)
# cell_4 = (0,0,0,1)
# cell_5 = (0,0,0,1)
# cell_6 = (0,0,0,1)
# cell_7 = (0,0,0,1)
sce <- SingleCellExperiment(list(rep1=cm))
noname_sce<- SingleCellExperiment(cm)
invalid_sce <- SingleCellExperiment(list(rep1=cm, rep2=cm))
perm_p_s = compute_permp(x=sce,
                         cluster_info=clusters,
                         perm.size=10,
                         bin_type="square",
                         bin_param=c(2,2),
                         test_genes=row.names(cm),
                         correlation_method = "pearson",
                         n_cores=1,
                         correction_method="BH",
                         w_x=w_x ,
                         w_y=w_y)
perm_noname= compute_permp(x=noname_sce,
                         cluster_info=clusters,
                         perm.size=10,
                         bin_type="square",
                         bin_param=c(2,2),
                         test_genes=row.names(cm),
                         correlation_method = "pearson",
                         n_cores=1,
                         correction_method="BH",
                         w_x=w_x ,
                         w_y=w_y)
test_that("Test permutation result - output dimension matches", {
    expect_equal(length(perm_p_s), 4)
    expect_equal(dim(perm_p_s$perm.arrays), c(4,1,10))
    expect_equal(dim(perm_p_s$obs.stat), c(4,1))
    expect_equal(dim(perm_p_s$perm.pval.adj), c(4,1))
    expect_equal(dim(perm_p_s$perm.pval), c(4,1))
    expect_equal(names(perm_hex_lst),
                 c("obs.stat", "perm.arrays", "perm.pval", "perm.pval.adj"))
    
})

test_that("Can work for one sample sce without sample name", {
    expect_equal(as.vector(perm_p_s$obs.stat),
                 as.vector(perm_noname$obs.stat))
    
})

test_that("Invaid input",{
    expect_error(compute_permp(x=invalid_sce,
                               cluster_info=clusters,
                               perm.size=10,
                               bin_type="square",
                               bin_param=c(2,2),
                               test_genes=row.names(cm),
                               correlation_method = "pearson",
                               n_cores=1,
                               correction_method="BH",
                               w_x=w_x ,
                               w_y=w_y))
})