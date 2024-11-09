library(SingleCellExperiment)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(SFEData)
library(Matrix)
library(BumpyMatrix)

##############################################################################
set.seed(200)
counts_sp1 <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
counts_sp2 <- matrix(rpois(100, lambda = 5), ncol=10, nrow=10)
##############################################################################

sce <- SingleCellExperiment(list(sp1=counts_sp1, sp2=counts_sp2))
sce_lst <- convert_data(sce, sample_names = c("sp1","sp2"))
exp_cm_lst <- sce@assays@data@listData
test_that("Invalid input",{
    expect_error(convert_data(sce, sample_names = c("sp01","sp2")))
})
test_that("Invalid input",{
    expect_error(convert_data("count", sample_names=c("sample1")))
    expect_error(convert_data(as.data.frame(counts_sp1),
                              sample_names=c("sample1")))
    expect_error(convert_data(counts_sp1, sample_names=c("sample1")))
    expect_error(convert_data(list("sample1" = counts_sp1),
                              sample_names=c("sample1")))
    expect_error(convert_data(list(counts_sp1),
                              sample_names=c("sample1")))
})

test_that("Test can create object from SingleCellExperiment", {
    expect_equal(sce_lst$trans_lst,NULL)
    expect_equal(length(sce_lst),2)
    expect_equal(names(sce_lst$cm_lst), c(c("sp1","sp2")))
    expect_equal(sce_lst$cm_lst, exp_cm_lst)
})
##############################################################################
# Visium barcode location from Space Ranger
data("visium_row_col")
coords1 <- visium_row_col[visium_row_col$col < 6 & visium_row_col$row < 6,]
coords1$row <- coords1$row * sqrt(3)

# Random toy sparse matrix
set.seed(29)
col_inds <- sample(1:13, 13)
row_inds <- sample(1:5, 13, replace = TRUE)
values <- sample(1:5, 13, replace = TRUE)
mat <- sparseMatrix(i = row_inds, j = col_inds, x = values)
colnames(mat) <- coords1$barcode
rownames(mat) <- sample(LETTERS, 5)
sfe3 <- SpatialFeatureExperiment(list(counts = mat), colData = coords1,
                                 spatialCoordsNames = c("col", "row"),
                                 spotDiameter = 0.7)
sfe3_example <- convert_data(sfe3, sample_names = c("sample01"))

test_that("Test can create object from SpatialFeatureExperiment", {
    expect_equal(sfe3_example$trans_lst,NULL)
    expect_equal(length(sfe3_example$cm_lst),1)
    expect_equal(names(sfe3_example$cm_lst), "sample01")
    expect_equal(as.vector(sfe3_example$cm_lst$sample01), as.vector(mat))
})

##############################################################################
# Construct toy dataset with 2 SFE samples
sfe1 <- McKellarMuscleData(dataset = "small")
sfe2 <- McKellarMuscleData(dataset = "small2")
spotPoly(sfe2)$sample_id <- "sample02"
sfe_combined <- SpatialFeatureExperiment::cbind(sfe1, sfe2)
sfe_lsts <- convert_data(sfe_combined, sample_names = c("Vis5A","sample02"))

test_that("Test can create object from multisample SFE", {
    expect_equal(sfe_lsts$trans_lst,NULL)
    expect_equal(length(sfe_lsts$cm_lst),2)
    expect_equal(names(sfe_lsts$cm_lst), c("Vis5A","sample02"))
})
##############################################################################
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
        molecules = mol))
test_that("Invalid input",{
    expect_error(convert_data(spe, sample_names="sample1"))
})
    
onesp_example <- convert_data(spe, sample_names="sample01")
prov_sum <- as.data.frame(table(df$gene))
cal_sum <- table(onesp_example$trans_lst$sample01$feature_name)
cal_sum <- as.data.frame(cal_sum)
prov_sum <- prov_sum[match(gs,prov_sum$Var1),]
cal_sum <- cal_sum[match(gs,cal_sum$Var1,),]

test_that("Can work for SE input",{
    expect_equal(names(onesp_example$trans_lst),"sample01")
    expect_equal(names(onesp_example$cm_lst),"sample01")
    expect_equal(onesp_example$cm_lst$sample01[,1],y[,1])
    expect_equal(onesp_example$cm_lst$sample01[1,],y[1,])
    expect_equal(prov_sum$Freq, cal_sum$Freq )
})

##############################################################################
# multi-sample SPE object
spe2 <- spe1 <-spe
spe1$sample_id <- paste(spe1$sample_id, "A", sep = ".")
spe2$sample_id <- paste(spe2$sample_id, "B", sep = ".")

# combine into single object
spe3 <- SpatialFeatureExperiment::cbind(spe1, spe2)

test_that("Invalid input",{
    expect_error(convert_data(spe, sample_names=c("sample1","sp2")))
})

