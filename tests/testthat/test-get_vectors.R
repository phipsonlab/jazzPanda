# simulate coordinates for clusters
clusters = data.frame(x = c(1,2,20,21,22,23,24),
                      y = c(23, 24, 1,2,3,4,5), cluster="A")
clusters$sample="rep1"

w_x=c(0,25)
w_y=c(0,25)
vecs_lst = get_vectors(data_lst= NULL, cluster_info = clusters,
                       bin_type = "square",
                       bin_param = c(2,2),
                       all_genes = NULL,
                       w_x = w_x, w_y=w_y)
test_that("Test can only vectorise clusters - output length matches", {
    expect_equal(length(vecs_lst), 1)
})


test_that("Test can only vectorise clusters - output vector matches", {
    expect_equal(as.vector(vecs_lst$cluster_mt), c(2,0,0,5))
})

#############################################################################
# simulate coordinates for genes
trans = data.frame(x = c(1,2,20,21,22,23,24),
                   y = c(23, 24, 1,2,3,4,5),
                   feature_name="A")
data=list(trans_info=trans)
w_x=c(0,25)
w_y=c(0,25)
vecs_lst_gene = get_vectors(data_lst= list("rep1"=data),
                            cluster_info = NULL,
                            bin_type = "square",
                            bin_param = c(2,2),
                            all_genes = c("A"),
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
clusters = data.frame(x = c(3, 5,11,21,2,23,19),
                      y = c(20, 24, 1,2,3,4,5), cluster="cluster_1")
clusters$sample="rep1"

data=list(trans_info=trans)
w_x=c(0,25)
w_y=c(0,25)
vecs_lst_gene = get_vectors(data_lst= list("rep1"= data),
                            cluster_info = clusters,
                            bin_type = "square",
                            bin_param = c(2,2),
                            all_genes = c("A","B","C"),
                            w_x = w_x, w_y=w_y)


test_that("Test can vectorise genes and clusters - output mathces", {
    expect_equal(as.vector(vecs_lst_gene$cluster_mt), c(2,0,2,3))
    expect_equal(as.vector(vecs_lst_gene$gene_mt), c(2, 0, 0 ,5 ,
                                                     1,0, 0, 1,
                                                     2, 0, 0, 5))
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

vecs_lst = get_vectors(data_lst= NULL, cluster_info = clusters,
                       cm_lst=list(rep1=cm),
                       bin_type = "square",
                       bin_param = c(2,2),
                       all_genes = row.names(cm),
                       w_x = w_x, w_y=w_y)

test_that("Test can use count matrix (square) - output vector matches", {
    expect_equal(as.vector(vecs_lst$gene_mt),
                 c(0,0,0,4,
                   8,0,0,31,
                   5,0,0,14,
                   1,0,0,5))
    
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
vecs_lst_cm = get_vectors(data_lst= NULL, cluster_info = clusters,
                          cm_lst=list(rep1=cm),
                          bin_type = "square",
                          bin_param = c(2,2),
                          all_genes = row.names(cm),
                          w_x = w_x, w_y=w_y)
vecs_lst_tr = get_vectors(data_lst= list("rep1" = list(trans_info = transcript_df)),
                          cluster_info = clusters,
                          bin_type = "square",
                          bin_param = c(2,2),
                          all_genes = row.names(cm),
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
vecs_lst = get_vectors(data_lst= NULL, cluster_info = clusters,
                       cm_lst=list(rep1=cm),
                       bin_type = "hexagon",
                       bin_param = c(10),
                       all_genes = row.names(cm),
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
clusters = data.frame(x = c(3, 5,11,21,2,23,19),
                      y = c(20, 24, 1,2,3,4,5), cluster="cluster_1")
clusters$sample="rep1"

data=list(trans_info=trans)
w_x=c(0,25)
w_y=c(0,25)
vecs_lst_gene = get_vectors(data_lst= list("rep1"= data),
                            cluster_info = clusters,
                            bin_type = "square",
                            bin_param = c(2,2),
                            all_genes = c("A","B","C"),
                            w_x = w_x, w_y=w_y,n_cores = 2)


test_that("Test can vectorise genes and clusters - output mathces", {
    expect_equal(as.vector(vecs_lst_gene$cluster_mt), c(2,0,2,3))
    expect_equal(as.vector(vecs_lst_gene$gene_mt), c(2, 0, 0 ,5 ,
                                                     1,0, 0, 1,
                                                     2, 0, 0, 5))
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
clusters = data.frame(x = c(3, 5,11,21,2,23,19),
                      y = c(20, 24, 1,2,3,4,5), cluster="cluster_1")
clusters$sample="rep1"

data=list(trans_info=trans)
w_x=c(0,25)
w_y=c(0,25)
vecs_lst_gene = get_vectors(data_lst= list("rep1"= data),
                            cluster_info = clusters,
                            bin_type = "square",
                            bin_param = c(2,2),
                            all_genes = c("A","B","C"),
                            w_x = w_x, w_y=w_y,n_cores = 2 )


test_that("Test can vectorise genes and clusters - output mathces", {
    expect_equal(as.vector(vecs_lst_gene$cluster_mt), c(2,0,2,3))
    expect_equal(as.vector(vecs_lst_gene$gene_mt), c(2, 0, 0 ,5 ,
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

head(coords_data)

vector_lst_1core = get_vectors(data_lst= list(sample1 = list(trans_info=coords_data)),
                                cluster_info = NULL,
                                bin_type="square",
                                bin_param=c(50,50),
                                all_genes = feature_names,
                                w_x=c(0,10000), w_y=c(0,10000), n_cores = 1)

vector_lst_5core = get_vectors(data_lst= list(sample1 =list(trans_info=coords_data)),
                               cluster_info = NULL,
                               bin_type="square",
                               bin_param=c(50,50),
                               all_genes = feature_names,
                               w_x=c(0,10000), w_y=c(0,10000), n_cores = 2)
test_that("Test can result from sequential matches with result from parallel", {
    expect_equal(as.vector(vector_lst_1core$gene_mt), 
                 as.vector(vector_lst_5core$gene_mt))

})
