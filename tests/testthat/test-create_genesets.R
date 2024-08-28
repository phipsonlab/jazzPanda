# simulate coordinates for genes
set.seed(15)

trans = as.data.frame(rbind(cbind(x = runif(10, min=1, max=10),
                                  y = runif(10, min=1, max=10), feature_name="A"),
                            cbind(x = runif(5, min=10, max=24),
                                  y = runif(5, min=1, max=10), feature_name="B"),
                            cbind(x = runif(10, min=10, max=24),
                                  y = runif(10, min=10, max=24), feature_name="C")))
trans$x = as.numeric(trans$x)
trans$y = as.numeric(trans$y)

data=list(trans_info=trans)
data_invalid = list(trans_info=trans, dummy_nme = list())
data_invalid_noname = list(trans)
invalid_trans_df = trans[,c("x","feature_name")]
data_invalid_missing_trans_cols = list(trans_info=invalid_trans_df)
geneset_res_A = create_genesets(data_lst=list("rep1"= data),
                                name_lst=list(dummy_W=c("A")),
                                bin_type="square",
                                bin_param=c(2,2),
                                w_x=c(0,25), w_y=c(0,25))
test_that("Invaid input",{
expect_error(create_genesets(data_lst=list("rep1"= data),
                             name_lst=list(dummy_W=c("A")),
                             bin_type="hexagon",
                             bin_param=c(2,2),
                             w_x=c(0,25), w_y=c(0,25)))
expect_error(create_genesets(data_lst=list("rep1"= data),
                             name_lst=list(dummy_W=c("A")),
                             bin_type="circle",
                             bin_param=c(2,2),
                             w_x=c(0,25), w_y=c(0,25)))
expect_error(create_genesets(data_lst=list("rep1"= data),
                             name_lst=list(dummy_W=c("A")),
                             bin_type="square",
                             bin_param=c(2),
                             w_x=c(0,25), w_y=c(0,25)))
# Line 47 
expect_error(create_genesets(data_lst=data_invalid,
                             name_lst=list(dummy_W=c("A")),
                             bin_type="square",
                             bin_param=c(2,2),
                             w_x=c(0,25), w_y=c(0,25)))

expect_error(create_genesets(data_lst=list("rep1"= data_invalid_noname),
                             name_lst=list(dummy_W=c("A")),
                             bin_type="square",
                             bin_param=c(2,2),
                             w_x=c(0,25), w_y=c(0,25)))

expect_error(create_genesets(data_lst=list("rep1"= data_invalid_missing_trans_cols),
                             name_lst=list(dummy_W=c("A")),
                             bin_type="square",
                             bin_param=c(2,2),
                             w_x=c(0,25), w_y=c(0,25)))

})
test_that("Test can create vectors for a single gene set with a single element A", {
    expect_equal(colnames(geneset_res_A), c("dummy_W"))
    expect_equal(as.vector(geneset_res_A$dummy_W), c(0,0,10, 0))
})
geneset_res_C = create_genesets(data_lst=list("rep1"= data),
                                name_lst=list(dummy_W=c("C")),
                                bin_type="square",
                                bin_param=c(2,2),
                                w_x=c(0,25), w_y=c(0,25))

test_that("Test can create vectors for a single gene set with a single element C", {
    expect_equal(colnames(geneset_res_C), c("dummy_W"))
    expect_equal(as.vector(geneset_res_C$dummy_W), c(2,7,0,1))
})

geneset_res1 = create_genesets(data_lst=list("rep1"= data),
                               name_lst=list(dummy_W=c("A","B")),
                               bin_type="square",
                               bin_param=c(2,2),
                               w_x=c(0,25), w_y=c(0,25))

test_that("Test can create vectors for a single gene set with multiple elements", {
    expect_equal(colnames(geneset_res1), c("dummy_W"))
    expect_equal(as.vector(geneset_res1$dummy_W), c(0,0,11,4))
})

geneset_res2 = create_genesets(data_lst=list("rep1"= data),
                               name_lst=list(dummy_A=c("A","C"),
                                             dummy_B=c("A","B","C")),
                               bin_type="square",
                               bin_param=c(2,2),
                               w_x=c(0,25), w_y=c(0,25))

test_that("Test can create vectors for gene sets- output mathces", {
    expect_equal(colnames(geneset_res2), c("dummy_A","dummy_B"))
    expect_equal(as.vector(geneset_res2$dummy_A), c(2,7,10,1))
    expect_equal(as.vector(geneset_res2$dummy_B), c(2,7,11,5))
})


geneset_res3 = create_genesets(data_lst=list("rep1"= data),
                               name_lst=list(dummy_A=c("A","B","C","C"),
                                             dummy_B=c("A","B","C")),
                               bin_type="square",
                               bin_param=c(2,2),
                               w_x=c(0,25), w_y=c(0,25))

test_that("Test can create vectors for duplicated gene sets", {
    expect_equal(colnames(geneset_res3), c("dummy_A","dummy_B"))
    expect_equal(as.vector(geneset_res3$dummy_A), c(2,7,11,5))
    expect_equal(as.vector(geneset_res3$dummy_B), c(2,7,11,5))
})

geneset_res4 = create_genesets(data_lst=list("rep1"= data),
                               name_lst=list(dummy_A=c("F","G"),
                                             dummy_B=c("A","B","C")),
                               bin_type="square",
                               bin_param=c(2,2),
                               w_x=c(0,25), w_y=c(0,25))

test_that("Test can create vectors for gene sets with non-overlapped genes", {
    expect_equal(colnames(geneset_res4), c("dummy_A","dummy_B"))
    expect_equal(as.vector(geneset_res4$dummy_A), c(0,0,0,0))
    expect_equal(as.vector(geneset_res4$dummy_B), c(2,7,11,5))
})

################
# define gene sets from count matrix
cm <- data.frame(rbind("gene_A"=c(0,0,2,0,0,0,2),
                       "gene_B"=c(5,3,3,13,0,1,14),
                       "gene_C"=c(5,0,1,5,1,0,7),
                       "gene_D"=c(0,1,1,2,0,0,2)))
colnames(cm)= paste("cell_", 1:7, sep="")

# simulate coordiantes for clusters
clusters = data.frame(x = c(1, 2,20,21,22,23,24),
                      y = c(23, 24, 1,2,3,4,5), cluster="A")
clusters$sample="sample1"
clusters$cell_id= colnames(cm)
invalid_cluster = clusters[,c("x","y")]
# simulate coordiantes for genes
w_x=c(0,25)
w_y=c(0,25)
geneset_res_cm = create_genesets(data_lst=list("sample1"= cm),
                               name_lst=list(dummy_A=c("gene_A","gene_B"),
                                             dummy_B=c("gene_C","gene_D")),
                               bin_type="square",
                               bin_param=c(2,2),
                               cluster_info = clusters,
                               w_x=w_x, w_y=w_y)

test_that("Test can define gene sets from count matrix", {
    expect_equal(dim(geneset_res_cm), c(4,2))
    expect_equal(as.vector(geneset_res_cm$dummy_A), c(8,0,0,35))
    expect_equal(as.vector(geneset_res_cm$dummy_B), c(6,0,0,19))
    expect_error(create_genesets(data_lst=list("sample1"= cm),
                                 name_lst=list(dummy_A=c("gene_A","gene_B"),
                                               dummy_B=c("gene_C","gene_D")),
                                 bin_type="square",
                                 bin_param=c(2,2),
                                 cluster_info = NULL,
                                 w_x=w_x, w_y=w_y))
    expect_error(create_genesets(data_lst=list("sample1"= cm),
                                 name_lst=list(dummy_A=c("gene_A","gene_B"),
                                               dummy_B=c("gene_C","gene_D")),
                                 bin_type="square",
                                 bin_param=c(2,2),
                                 cluster_info = invalid_cluster,
                                 w_x=w_x, w_y=w_y))
})

### check input
# simulate coordinates for genes
set.seed(15)
trans = as.data.frame(rbind(cbind(x = runif(10, min=1, max=10),
                                  y = runif(10, min=1, max=10), feature_name="A"),
                            cbind(x = runif(5, min=10, max=24),
                                  y = runif(5, min=1, max=10), feature_name="B"),
                            cbind(x = runif(10, min=10, max=24),
                                  y = runif(10, min=10, max=24), feature_name="C")))
trans$x = as.numeric(trans$x)
trans$y = as.numeric(trans$y)

data=list(trans_info=trans)

check_res1 = check_geneset_input(data_lst=list("rep1"= data),
                                 bin_type="square",
                                 bin_param=c(2,2),
                                 w_x=c(0,25), w_y=c(0,25),cluster_info=NULL)

test_that("Test can validate the correct input - square bins", {
    expect_equal(check_res1$bin_length, 4)
    expect_equal(check_res1$use_cm,FALSE)
})


check_res2 = check_geneset_input(data_lst=list("rep1"= data),
                                 bin_type="hexagon",
                                 bin_param=c(2),
                                 w_x=c(0,25), w_y=c(0,25),
                                 cluster_info=NULL)

test_that("Test can validate the correct input - hex bins", {
    expect_equal(check_res2$bin_length, 77)
    expect_equal(check_res2$use_cm,FALSE)
})


check_res3 = check_geneset_input(data_lst=list("rep1"= data),
                                 bin_type="hexagon",
                                 bin_param=c(5),
                                 w_x=c(0,25), w_y=c(0,25),
                                 cluster_info=NULL)

test_that("Test can validate the correct input - hex bins", {
    expect_equal(check_res3$bin_length, 17)
    expect_equal(check_res3$use_cm, FALSE)
})
