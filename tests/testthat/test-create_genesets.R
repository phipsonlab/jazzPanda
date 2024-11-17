library(SpatialExperiment)
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
trans$cell =  sample(paste("cell",1:7, sep=""), 
                     replace = TRUE, size = nrow(trans))
mol <- BumpyMatrix::splitAsBumpyMatrix(
    trans[, c("x", "y")], 
    row = trans$feature_name, col = trans$cell )

spe_rep1 <- SpatialExperiment(
    assays = list(molecules = mol),sample_id ="rep1" )

data=trans

invalid_trans_df = trans[,c("x","feature_name","cell")]
mol_inv <- BumpyMatrix::splitAsBumpyMatrix(
    invalid_trans_df[, c("x")], 
    row = invalid_trans_df$feature_name, col = invalid_trans_df$cell )

spe_invalid <- SpatialExperiment(
    assays = list(molecules = mol_inv),sample_id ="rep1" )


#data_invalid_missing_trans_cols = list(trans_info=invalid_trans_df)
geneset_res_A = create_genesets(x=spe_rep1,sample_names = "rep1",
                                name_lst=list(dummy_W=c("A")),
                                bin_type="square",
                                bin_param=c(2,2),cluster_info = NULL,
                                w_x=c(0,25), w_y=c(0,25))
geneset_res_A_hex = create_genesets(x=spe_rep1,sample_names = "rep1",
                                    name_lst=list(dummy_W=c("A")),
                                bin_type="hexagon",
                                bin_param=c(10),cluster_info = NULL,
                                w_x=c(0,25), w_y=c(0,25))
test_that("Invaid input",{
expect_error(create_genesets(x=spe_rep1,sample_names = "rep1",
                             name_lst=list(dummy_W=c("A")),
                             bin_type="hexagon",
                             bin_param=c(2,2),cluster_info = NULL,
                             w_x=c(0,25), w_y=c(0,25)))
expect_error(create_genesets(x=spe_rep1,sample_names = "rep1",
                             name_lst=list(dummy_W=c("A")),
                             bin_type="circle",
                             bin_param=c(2,2),cluster_info = NULL,
                             w_x=c(0,25), w_y=c(0,25)))
expect_error(create_genesets(x=spe_rep1,sample_names = "rep1",
                             name_lst=list(dummy_W=c("A")),
                             bin_type="square",
                             bin_param=c(2),cluster_info = NULL,
                             w_x=c(0,25), w_y=c(0,25)))

})
test_that("Test can create vectors for a single gene set with a single element A", {
    expect_equal(colnames(geneset_res_A), c("dummy_W"))
    expect_equal(as.vector(geneset_res_A$dummy_W), c(0,0,10, 0))
    expect_equal(as.vector(geneset_res_A_hex$dummy_W), c(2,5,0,3,0,0,0))
})

geneset_res_C = create_genesets(x=spe_rep1,sample_names = "rep1",
                                name_lst=list(dummy_W=c("C")),
                                bin_type="square",
                                bin_param=c(2,2),cluster_info = NULL,
                                w_x=c(0,25), w_y=c(0,25))

test_that("Test can create vectors for a single gene set with a single element C", {
    expect_equal(colnames(geneset_res_C), c("dummy_W"))
    expect_equal(as.vector(geneset_res_C$dummy_W), c(2,7,0,1))
})

geneset_res1 = create_genesets(x=spe_rep1,sample_names = "rep1",
                               name_lst=list(dummy_W=c("A","B")),
                               bin_type="square",
                               bin_param=c(2,2),cluster_info = NULL,
                               w_x=c(0,25), w_y=c(0,25))

test_that("Test can create vectors for a single gene set with multiple elements", {
    expect_equal(colnames(geneset_res1), c("dummy_W"))
    expect_equal(as.vector(geneset_res1$dummy_W), c(0,0,11,4))
})

geneset_res2 = create_genesets(x=spe_rep1,sample_names = "rep1",
                               name_lst=list(dummy_A=c("A","C"),
                                             dummy_B=c("A","B","C")),
                               bin_type="square",
                               bin_param=c(2,2),cluster_info = NULL,
                               w_x=c(0,25), w_y=c(0,25))

test_that("Test can create vectors for gene sets- output mathces", {
    expect_equal(colnames(geneset_res2), c("dummy_A","dummy_B"))
    expect_equal(as.vector(geneset_res2$dummy_A), c(2,7,10,1))
    expect_equal(as.vector(geneset_res2$dummy_B), c(2,7,11,5))
})


geneset_res3 = create_genesets(x=spe_rep1,sample_names = "rep1",
                               name_lst=list(dummy_A=c("A","B","C","C"),
                                             dummy_B=c("A","B","C")),
                               bin_type="square",cluster_info = NULL,
                               bin_param=c(2,2),
                               w_x=c(0,25), w_y=c(0,25))

test_that("Test can create vectors for duplicated gene sets", {
    expect_equal(colnames(geneset_res3), c("dummy_A","dummy_B"))
    expect_equal(as.vector(geneset_res3$dummy_A), c(2,7,11,5))
    expect_equal(as.vector(geneset_res3$dummy_B), c(2,7,11,5))
})

geneset_res4 = create_genesets(x=spe_rep1,sample_names = "rep1",
                               name_lst=list(dummy_B=c("A","B","C")),
                               bin_type="square",
                               bin_param=c(2,2),cluster_info = NULL,
                               w_x=c(0,25), w_y=c(0,25))

test_that("Test can not create vectors for gene sets with non-overlapped genes", {
    expect_error(create_genesets(x=spe_rep1,sample_names = "rep1",
                                 name_lst=list(dummy_A=c("F","G")),
                                 bin_type="square",
                                 bin_param=c(2,2),cluster_info = NULL,
                                 w_x=c(0,25), w_y=c(0,25)))
    expect_equal(colnames(geneset_res4), c("dummy_B"))
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
sce <- SingleCellExperiment(list(sample1=cm))
# simulate coordiantes for genes
w_x=c(0,25)
w_y=c(0,25)
geneset_res_cm = create_genesets(x=sce,sample_names = "sample1",
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
    expect_error(create_genesets(x=sce,sample_names = "sample1",
                                 name_lst=list(dummy_A=c("gene_A","gene_B"),
                                               dummy_B=c("gene_C","gene_D")),
                                 bin_type="square",
                                 bin_param=c(2,2),
                                 cluster_info = invalid_cluster,
                                 w_x=w_x, w_y=w_y))
})

############################################################################
# example of two SPE objects but with different gene sets 
# real world exmple: if the two samples have slightly different negative controls probes
s1 =as.data.frame(cbind(feature_name = c("A","B","C"),
           x= c(1,2,3),
           y= c(24, 22, 23), 
           cell_id = c("cell1","cell2","cell3")))
s2 = as.data.frame(cbind(feature_name = c("A","E","F","G"),
           x= c(1,2,3,4),
           y= c(24, 22, 23,24), 
           cell_id = c("cell4","cell5","cell6","cell7")))

shared_nc =  intersect(s1$feature_name,s2$feature_name)

s1_added = cbind(x=NA, y=NA, 
                   feature_name=setdiff(s2$feature_name,shared_nc),
                   cell_id = NA)
s2_added = cbind(x=NA, y=NA, 
                   feature_name=setdiff(s1$feature_name,shared_nc),
                   cell_id = NA)

s1 = rbind(s1,s1_added)
s2 = rbind(s2,s2_added)

s1$category = "probe"
s1[s1$feature_name %in% c("A","E"), "category"] = "codeword"

s2$category = "probe"
s2[s2$feature_name %in% c("C","G"), "category"] = "codeword"

s1_mol <- BumpyMatrix::splitAsBumpyMatrix(
    s1[, c("x", "y","category")], 
    row = s1$feature_name, col = s1$cell_id )

s1_spe<- SpatialExperiment(
    assays = list(molecules = s1_mol),sample_id ="sample1")

s2_mol <- BumpyMatrix::splitAsBumpyMatrix(
    s2[, c("x", "y","category")], 
    row = s2$feature_name, col = s2$cell_id )

s2_spe<- SpatialExperiment(
    assays = list(molecules = s2_mol),sample_id ="sample2")

combined_spe = BiocGenerics::cbind(s1_spe, s2_spe)

vecs_lst = get_vectors(x=combined_spe,sample_names=c("sample1","sample2"), 
            bin_type="square",test_genes = c("A","B","C","E","F","G"),
            bin_param=c(2,2), 
            w_x=w_x, w_y=w_y,
            cluster_info = NULL)
test_that("Test can not create vectors for multi-samples with different genes", {
    expect_equal(as.vector(vecs_lst$gene_mt[,"A"]), c(1,0,0,0,1,0,0,0))
    expect_equal(as.vector(vecs_lst$gene_mt[,"B"]), c(1,0,0,0,0,0,0,0))
    expect_equal(as.vector(vecs_lst$gene_mt[,"C"]), c(1,0,0,0,0,0,0,0))
    expect_equal(as.vector(vecs_lst$gene_mt[,"E"]), c(0,0,0,0,1,0,0,0))
    expect_equal(as.vector(vecs_lst$gene_mt[,"F"]), c(0,0,0,0,1,0,0,0))
    expect_equal(as.vector(vecs_lst$gene_mt[,"G"]), c(0,0,0,0,1,0,0,0))
})


example_vec = create_genesets(x=combined_spe,
                                     sample_names=c("sample1","sample2"),
                                     name_lst=list(probe=c("A","E"), 
                                                   codeword=c("C","G")),
                                     bin_type="square",
                                     bin_param=c(2,2), 
                                     w_x=w_x, w_y=w_y,
                                     cluster_info = NULL)
test_that("Test can not create vectors for multi-samples with different genes", {
    expect_equal(as.vector(example_vec[,"probe"]), c(1,0,0,0,2,0,0,0))
    expect_equal(as.vector(example_vec[,"codeword"]), c(1,0,0,0,1,0,0,0))
})
