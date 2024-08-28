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

geneset_res_A = create_genesets(data_lst=list("rep1"= data),
                                name_lst=list(dummy_W=c("A")),
                                bin_type="square",
                                bin_param=c(2,2),
                                w_x=c(0,25), w_y=c(0,25))

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
