library(SingleCellExperiment)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(SFEData)
library(Matrix)
library(BumpyMatrix)

# simulate coordinates for clusters
clusters = data.frame(x = c(1,2,20,21,22,23,24),
                      y = c(23, 24, 1,2,3,4,5), cluster="A")
clusters$sample="rep1"

w_x=c(0,25)
w_y=c(0,25)
vecs_lst = get_vectors(x= NULL, cluster_info = clusters,
                       bin_type = "square",sample_names="rep1",
                       bin_param = c(2,2),
                       test_genes = NULL,
                       w_x = w_x, w_y=w_y)

test_that("Test can only vectorise clusters - output length matches", {
    expect_equal(length(vecs_lst), 1)
})


test_that("Test can only vectorise clusters - output vector matches", {
    expect_equal(as.vector(vecs_lst$cluster_mt), c(2,0,0,5))
})

#############################################################################
# one sample SPE object
n <- 1e3  # number of molecules
ng <- 50  # number of genes
nc <- 20  # number of cells
# sample xy-coordinates in [0, 1]
x <- runif(n)
y <- runif(n)
# assign each molecule to some gene-cell pair
gs <- paste0("gene", seq(ng))
cs <- paste0("cell", seq(nc))
gene <- sample(gs, n, TRUE)
cell <- sample(cs, n, TRUE)
# assure gene & cell are factors so that
# missing observations aren't dropped
gene <- factor(gene, gs)
cell <- factor(cell, cs)
# construct data.frame of molecule coordinates
df <- data.frame(gene, cell, x, y)
mol <- BumpyMatrix::splitAsBumpyMatrix(
    df[, c("x", "y")], 
    row = gene, col = cell)
y <- with(df, table(gene, cell))
y <- as.matrix(unclass(y))
# construct SpatialExperiment
spe <- SpatialExperiment(
    assays = list(
        counts = y, 
        molecules = mol),sample_id ="sample1" )
vecs_lst_gene = get_vectors(x= spe,
                            sample_names = c("sample1"),
                            cluster_info = NULL,
                            bin_type = "square",
                            bin_param = c(2,2),
                            test_genes = row.names(spe),
                            w_x = w_x, w_y=w_y)
test_that("Test can only vectorise genes - output length mathces", {
    expect_equal(length(vecs_lst_gene), 1)
})

#############################################################################
# simulate coordinates for genes
trans = as.data.frame(rbind(cbind(x = c(1,2,20,21,22,23,24),
                                  y = c(23, 24, 1,2,3,4,5),
                                  feature_name="A"),
                            cbind(x = c(1,20),
                                  y = c(15, 10),
                                  feature_name="B"),
                            cbind(x = c(1,2,20,21,22,23,24),
                                  y = c(23, 24, 1,2,3,4,5),
                                  feature_name="C")))

trans$x = as.numeric(trans$x)
trans$y = as.numeric(trans$y)
trans$cell = sample(c("cell1", "cell2","cell3"), 
                    replace = TRUE, size = nrow(trans))
mol <- BumpyMatrix::splitAsBumpyMatrix(
    trans[, c("x", "y")], 
    row = trans$feature_name, col = trans$cell )

clusters = data.frame(x = c(3, 5,11,21,2,23,19),
                      y = c(20, 24, 1,2,3,4,5), cluster="cluster_1")
clusters$sample="sample1"
invalid_cluster = clusters[, c("x","y")]
invalid_data=trans[, c("x","y")]
invalid_trans_df = trans[,c("x","feature_name","cell")]
mol_inv <- BumpyMatrix::splitAsBumpyMatrix(
    invalid_trans_df[, c("x")], 
    row = invalid_trans_df$feature_name, col = invalid_trans_df$cell )

spe_invalid <- SpatialExperiment(
    assays = list(molecules = mol_inv),sample_id ="rep1" )


w_x=c(0,25)
w_y=c(0,25)
vecs_lst_cluster = get_vectors(x= NULL,
                            cluster_info = clusters,
                            sample_names="sample1",
                            bin_type = "square",
                            bin_param = c(2,2),
                            test_genes = NULL,
                            w_x = w_x, w_y=w_y)
spe <- SpatialExperiment(
    assays = list(molecules = mol),sample_id ="sample1" )

vecs_lst_gene = get_vectors(x= spe,sample_names="sample1" ,
                               cluster_info = NULL,
                               bin_type = "square",
                               bin_param = c(2,2),
                               test_genes = row.names(spe),
                               w_x = w_x, w_y=w_y)

test_that("Test can vectorise genes and clusters - output mathces", {
    expect_equal(as.vector(vecs_lst_cluster$cluster_mt), c(2,0,2,3))
    expect_equal(as.vector(vecs_lst_gene$gene_mt), c(2, 0, 0 ,5 ,
                                                     1,0, 0, 1,
                                                     2, 0, 0, 5))
})
test_that("Invaid input",{
    expect_error(get_vectors(x=spe,  sample_names="sample1",
                             cluster_info = clusters,
                             bin_type="square",
                             bin_param=c(2,2),
                             test_genes =matrix(0, nrow=5, ncol=2),
                             w_x=c(0,25), w_y=c(0,25)))
    
    expect_error(get_vectors(x=trans,  sample_names="sample1",
                             cluster_info = clusters,
                             bin_type="hexagon",
                             bin_param=c(2,2),
                             test_genes = row.names(spe),
                             w_x=c(0,25), w_y=c(0,25)))
    expect_error(get_vectors(x=NULL,  sample_names="sample1",
                                cluster_info = clusters,
                                bin_type="hexagon",
                                bin_param=c(2,2),
                                test_genes = row.names(spe),
                                w_x=c(0,25), w_y=c(0,25)))
    expect_error(get_vectors(x=NULL,  sample_names="sample1",
                                cluster_info = clusters,
                                bin_type="circle",
                                bin_param=c(2,2),
                                test_genes = row.names(spe),
                                w_x=c(0,25), w_y=c(0,25)))
    expect_error(get_vectors(x=NULL,  sample_names="sample1",
                                cluster_info = clusters,
                                bin_type="square",
                                bin_param=c(2),
                                w_x=c(0,25), w_y=c(0,25)))
    expect_error(get_vectors(x=NULL, cluster_info = NULL,
                             sample_names="sample1",
                             bin_type="hexagon",
                             bin_param=c(2),
                             test_genes = row.names(spe),
                             w_x=c(0,25), w_y=c(0,25)))
    expect_error(get_vectors(x=NULL,  sample_names="sample1",
                            cluster_info=invalid_cluster,
                            bin_type="square",
                            bin_param=c(2,2),
                            test_genes = row.names(spe),
                            w_x=c(0,25), w_y=c(0,25)))
    expect_error(get_vectors(x=spe,  sample_names="sample1",
                             cluster_info=invalid_cluster,
                             bin_type="square",
                             bin_param=c(2,2),
                             test_genes = row.names(spe),
                             w_x=c(0,25), w_y=c(0,25), use_cm=TRUE))
    
    #  can not match test_genes from x
    expect_error(get_vectors(x= spe,
                             sample_names="sample1",
                             test_genes =c("geneA","geneB"),
                             cluster_info=clusters,
                             bin_type="square",
                             bin_param=c(2,2),
                             w_x=c(0,25), w_y=c(0,25)))
    # mismatch sample_names and clsuetr_info$sample
    expect_error(get_vectors(x= NULL,
                             cluster_info = clusters,
                             sample_names="sample2",
                             bin_type = "square",
                             bin_param = c(2,2),
                             test_genes = NULL,
                             w_x = w_x, w_y=w_y))
    # missing y coordiantes for transcript 
    expect_error(get_vectors(x= spe_invalid,
                             cluster_info = NULL,
                             sample_names="rep1",
                             bin_type = "square",
                             bin_param = c(2,2),
                             test_genes = unique(invalid_trans_df$feature_name),
                             w_x = w_x, w_y=w_y))
    
})

#############################################################################
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
vecs_lst = get_vectors(x= sce, cluster_info = clusters,
                       sample_names = "rep1",
                       bin_type = "square",
                       bin_param = c(2,2),
                       test_genes = row.names(cm),
                       w_x = w_x, w_y=w_y)

test_that("Test can use count matrix (square) - output vector matches", {
    expect_equal(as.vector(vecs_lst$gene_mt),
                 c(0,0,0,4,
                   8,0,0,31,
                   5,0,0,14,
                   1,0,0,5))
    
})
#  Missing cluster information to build gene vector matrix.
test_that("Invaid input",{
    expect_error(get_vectors(x=sce,sample_names = "rep1",
                                cluster_info = NULL,
                                bin_type="square",
                                bin_param=c(2,2),test_genes = row.names(cm),
                                w_x=c(0,25), w_y=c(0,25)))
    
})

#############################################################################
# gene vector from count matrix matches with vector from transcript coordinates
cm <- data.frame(rbind("gene_A"=c(0,0,2,0,0,0,2),
                       "gene_B"=c(5,3,3,13,0,1,14),
                       "gene_C"=c(5,0,1,5,1,0,7),
                       "gene_D"=c(0,1,1,2,0,0,2)))
colnames(cm)= paste("cell_", 1:7, sep="")

# simulate coordinates for clusters
clusters = data.frame(x = c(1, 2,20,21,22,23,24),
                      y = c(23, 24, 1,2,3,4,5), cluster="A")
clusters$sample="rep1"
clusters$cell_id= colnames(cm)
# simulate coordinates for genes
w_x=c(0,25)
w_y=c(0,25)
# cell_1 = (1,0,0,0)
# cell_2 = (1,0,0,0)
# cell_3 = (0,0,0,1)
# cell_4 = (0,0,0,1)
# cell_5 = (0,0,0,1)
# cell_6 = (0,0,0,1)
# cell_7 = (0,0,0,1)
# if contains duplicated x,y pairs in transcript detections, then will be removed 
#       --> will show difference compared to count matrix 
# set.seed(12)
# transcript_df = as.data.frame(rbind(cbind(x = sample(12:24, size = 4),
#                                           y = sample(1:11, size = 4),
#                                           feature_name="gene_A"),
#                                     cbind(x = c(sample(1:11, size = 8,replace = TRUE),sample(12:24, size = 31,replace = TRUE)),
#                                           y = c(sample(14:24, size = 8,replace = TRUE),sample(1:11, size = 31,replace = TRUE)),
#                                           feature_name="gene_B"),
#                                     cbind(x = c(sample(1:11, size = 5),sample(14:24, size = 14,replace = TRUE)),
#                                           y = c(sample(14:24, size = 5),sample(1:11, size = 14,replace = TRUE)),
#                                           feature_name="gene_C"),
#                                     cbind(x = c(sample(1:11, size = 1),sample(14:24, size = 5,replace = TRUE)),
#                                           y = c(sample(14:24, size = 1),sample(1:11, size = 5,replace = TRUE)),
#                                           feature_name="gene_D")
# ))
set.seed(12)
transcript_df = as.data.frame(rbind(cbind(x = runif(min=12, max=24, n = 4),
                                          y = runif(min=1, max=11, n = 4),
                                          feature_name="gene_A"),
                                    cbind(x = c(runif(min=1, max=11,n = 8),runif(min=12, max=24, n = 31)),
                                          y = c(runif(min=14, max=24,n = 8),runif(min=1, max=11, n = 31)),
                                          feature_name="gene_B"),
                                    cbind(x = c(runif(min=1, max=11,n = 5),runif(min=14, max=24, n = 14)),
                                          y = c(runif(min=14, max=24, n = 5),runif(min=1, max=11,n = 14)),
                                          feature_name="gene_C"),
                                    cbind(x = c(runif(min=1, max=11,n = 1),runif(min=14, max=24, n = 5)),
                                          y = c(runif(min=14, max=24, n = 1),runif(min=1, max=11, n = 5)),
                                          feature_name="gene_D")
))

transcript_df$x=as.numeric(transcript_df$x)
transcript_df$y=as.numeric(transcript_df$y)
sce_m <- SingleCellExperiment(list(rep1=cm))
transcript_df$cell =  sample(paste("cell",1:7, sep=""), 
                             replace = TRUE, size = nrow(transcript_df))
mol <- BumpyMatrix::splitAsBumpyMatrix(
    transcript_df[, c("x", "y")], 
    row = transcript_df$feature_name, col = transcript_df$cell )

spe_m <- SpatialExperiment(
    assays = list(molecules = mol),sample_id ="rep1" )

vecs_lst_cm = get_vectors(x= sce_m, cluster_info = clusters,
                          sample_names = "rep1",
                          bin_type = "square",
                          bin_param = c(2,2),
                          test_genes = row.names(cm),
                          w_x = w_x, w_y=w_y)
vecs_lst_tr = get_vectors(x= spe_m,
                          cluster_info = clusters,
                          sample_names = "rep1",
                          bin_type = "square",
                          bin_param = c(2,2),
                          test_genes = row.names(cm),
                          w_x = w_x, w_y=w_y)
test_that("Test count matrix vectors match with transcript vectors - gene vectors", {
    expect_equal(as.vector(vecs_lst_tr$gene_mt),
                 as.vector(vecs_lst_cm$gene_mt))
    
})
test_that("Test count matrix vectors match with transcript vectors - cluster vectors", {
    expect_equal(as.vector(vecs_lst_tr$cluster_mt),
                 as.vector(vecs_lst_cm$cluster_mt))
    
})
# cell_1 = c()
vecs_lst = get_vectors(x= sce_m, cluster_info = clusters,
                       sample_names = "rep1",
                       bin_type = "hexagon",
                       bin_param = c(10),
                       test_genes = row.names(cm),
                       w_x = w_x, w_y=w_y)

test_that("Test can use count matrix (hex) - output vector matches", {
    expect_equal(as.vector(vecs_lst$gene_mt),
                 c(0,0,0,0,4,0,0,
                   0,0,0,0,31,8,0,
                   0,0,0,0,14,5,0,
                   0,0,0,0,5,1,0))
    
})
#############################################################################
# simulate coordiantes for genes with 3 cores 
trans = as.data.frame(rbind(cbind(x = c(1,2,20,21,22,23,24),
                                  y = c(23, 24, 1,2,3,4,5),
                                  feature_name="A"),
                            cbind(x = c(1,20),
                                  y = c(15, 10),
                                  feature_name="B"),
                            cbind(x = c(1,2,20,21,22,23,24),
                                  y = c(23, 24, 1,2,3,4,5),
                                  feature_name="C")))

trans$x = as.numeric(trans$x)
trans$y = as.numeric(trans$y)
trans$cell =  sample(paste("cell",1:7, sep=""), 
                             replace = TRUE, size = nrow(trans))
mol <- BumpyMatrix::splitAsBumpyMatrix(
    trans[, c("x", "y")], 
    row = trans$feature_name, col = trans$cell )

spe_rep1 <- SpatialExperiment(
    assays = list(molecules = mol),sample_id ="rep1" )


clusters = data.frame(x = c(3, 5,11,21,2,23,19),
                      y = c(20, 24, 1,2,3,4,5), cluster="cluster_1")
clusters$sample="rep1"
w_x=c(0,25)
w_y=c(0,25)
vecs_lst_full = get_vectors(x=spe_rep1,sample_names = "rep1",
                            cluster_info = clusters,
                            bin_type = "square",
                            bin_param = c(2,2),
                            test_genes = c("A","B","C"),
                            w_x = w_x, w_y=w_y,n_cores = 2)


test_that("Test can vectorise genes and clusters - output mathces", {
    expect_equal(as.vector(vecs_lst_full$cluster_mt), c(2,0,2,3))
    expect_equal(as.vector(vecs_lst_full$gene_mt), c(2, 0, 0 ,5 ,
                                                     1,0, 0, 1,
                                                     2, 0, 0, 5))
})

#############################################################################
# simulate coordiantes for genes with 3 cores 
# Set seed for reproducibility
set.seed(123)

# Generate 100 feature names
feature_names <- paste("gene", 1:200, sep="")

# Function to generate random coordinates for a given gene
generate_coords <- function(feature_name) {
    # Random number of coordinates between 100 and 10000
    n_coords <- sample(5000:10000, 1)
    
    # Generate random x and y coordinates
    x_coords <- runif(n_coords, min=0, max=10000)
    y_coords <- runif(n_coords, min=0, max=10000)
    
    # Create a data frame
    data.frame(x = x_coords, y = y_coords, feature_name = rep(feature_name, n_coords))
}

# Apply the function to each feature name and combine results
coords_data <- do.call(rbind, lapply(feature_names, generate_coords))

coords_data$cell =  sample(paste("cell",1:10, sep=""), 
                     replace = TRUE, size = nrow(coords_data))
mol <- BumpyMatrix::splitAsBumpyMatrix(
    coords_data[, c("x", "y")], 
    row = coords_data$feature_name, col = coords_data$cell )

spe_rep2 <- SpatialExperiment(
    assays = list(molecules = mol),sample_id ="rep2" )


vector_lst_1core = get_vectors(x= spe_rep2, sample_names = "rep2",
                                cluster_info = NULL,
                                bin_type="square",
                                bin_param=c(50,50),
                                test_genes = feature_names,
                                w_x=c(0,10000), w_y=c(0,10000), n_cores = 1)

vector_lst_5core = get_vectors(x= spe_rep2, sample_names = "rep2",
                               cluster_info = NULL,
                               bin_type="square",
                               bin_param=c(50,50),
                               test_genes = feature_names,
                               w_x=c(0,10000), w_y=c(0,10000), n_cores = 2)
test_that("Test can result from sequential matches with result from parallel",{
    expect_equal(as.vector(vector_lst_1core$gene_mt), 
                 as.vector(vector_lst_5core$gene_mt))

})
##############
# two sample SingleCellExperiment scenario
set.seed(200)
counts_sp1 <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
counts_sp2 <- matrix(rpois(100, lambda = 5), ncol=10, nrow=10)
colnames(counts_sp1) <- paste("cell", 1:10, sep="")
colnames(counts_sp2) <- paste("cell", 1:10, sep="")
row.names(counts_sp2) <- paste("gene", 1:10, sep="")
row.names(counts_sp1) <- paste("gene", 1:10, sep="")
clusters = rbind(data.frame(x = c(1, 2,20.5,21.8,22.6,23.5,24.4,23.3,24.2,24.1),
                      y = c(23.1, 24.2, 1,2,3,4,5,1,2,3)),
                 data.frame(x = c(1, 2,20,21,22,23,24,23,24,24),
                            y = c(23, 24, 1,2,3,4,5,1,2,3)))
clusters$sample=rep(c("sp1","sp2"), each=10)
clusters$cell_id= rep(paste("cell", 1:10, sep=""), time=2)
clusters$cluster = rep(c("A","B"), each=10)

# cell_1 = (1,0,0,0)
# cell_2 = (1,0,0,0)
# cell_3 = (0,0,0,1)
# cell_4 = (0,0,0,1)
# cell_5 = (0,0,0,1)
# cell_6 = (0,0,0,1)
# cell_7 = (0,0,0,1)
# cell_8 = (0,0,0,1)
# cell_9 = (0,0,0,1)
# cell_10 = (0,0,0,1)
sce_two <- SingleCellExperiment(list(sp1=counts_sp1, sp2=counts_sp2))


vector_lst_twosample = get_vectors(x= sce_two, sample_names =c("sp1","sp2"),
                               cluster_info = clusters,
                               bin_type="square",
                               bin_param=c(2,2),
                               test_genes = paste("gene", 1:10, sep=""),
                               w_x=c(0,25), w_y=c(0,25), 
                               n_cores = 1)
test_that("Test can result from sequential matches with result from parallel",{
    expect_equal(as.vector(vector_lst_twosample$cluster_mt[,"A"]), 
                 c(2,0,0,8,0,0,0,0))
    expect_equal(as.vector(vector_lst_twosample$cluster_mt[,"B"]), 
                 c(0,0,0,0,2,0,0,8))
    expect_equal(as.vector(vector_lst_twosample$gene_mt[,"gene1"]), 
                 c(17, 0,0,82, 12,0,0,42))
    
})


