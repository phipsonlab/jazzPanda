# Code to prepare data used in vignette

# The raw data can be accessed here:
# <https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast>
# <https://www.biorxiv.org/content/10.1101/2022.10.06.510405v1>
#  You can download the data with the following command:
#
# **Download input/output Files**
# - Replicate 1:
# wget <https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_gene_panel.json>
# wget <https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_panel.tsv>
# wget <https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_gene_groups.csv>
# wget <https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip>
# unzip Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip
# ```
# - Replicate 2:
# wget <https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep2/Xenium_FFPE_Human_Breast_Cancer_Rep2_gene_panel.json>
# wget <https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep2/Xenium_FFPE_Human_Breast_Cancer_Rep2_panel.tsv>
# wget <https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep2/Xenium_FFPE_Human_Breast_Cancer_Rep2_gene_groups.csv>
# wget <https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep2/Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip>
# unzip Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip
# ```
#
# Details about each file can be found here:
# <https://www.10xgenomics.com/support/in-situ-gene-expression/documentation/steps/onboard-analysis/understanding-xenium-outputs#overview>

# The data used in this vignette is the a small section of the two
# Xenium human breast cancer replicates.


library(data.table)
library(Seurat)
library(SpatialExperiment)
# download the official cell types of both samples
s1_url <- "https://raw.githubusercontent.com/phipsonlab/jazzPanda_vignette/main/Xenium_rep1_supervised_celltype.csv"
s2_url <- "https://raw.githubusercontent.com/phipsonlab/jazzPanda_vignette/main/Xenium_rep2_supervised_celltype.csv"
# Download the file
download.file(s1_url, "graphclust_rep1.csv")
download.file(s2_url, "graphclust_rep2.csv")

get_xenium_data<-function(path,mtx_name, trans_name="transcript_info.csv.gz",
                   cells_name="cell_info.csv.gz"){
    
    transcript_info <- as.data.frame(fread(paste(path, trans_name,sep="")))
    cell_info <- as.data.frame(fread(paste(path,cells_name,sep="")))
    
    data <- Read10X(data.dir = paste(path,mtx_name, sep=""))
    
    cm <- as.matrix(data$`Gene Expression`)
    r_codeword <- as.matrix(data$`Negative Control Codeword`)
    
    r_probe <- as.matrix(data$`Negative Control Probe`)
    cm_neg <- as.data.frame(rbind(r_probe, r_codeword))
    zero_cells <- colnames(cm)[colSums(cm)==0]
    
    transcript_info$x <- as.numeric(transcript_info$x_location)
    transcript_info$y <- as.numeric(transcript_info$y_location)
    transcript_info <- transcript_info[transcript_info$qv >=20 &
                                       transcript_info$cell_id != -1 &
                                       !(transcript_info$cell_id %in% zero_cells), ]
    cell_info$x <- as.numeric(cell_info$x_centroid)
    cell_info$y <- as.numeric(cell_info$y_centroid)
    
    return (list(trans_info=transcript_info[, c("x","y","feature_name","cell_id")], 
                 cell_info=cell_info, cm=cm,                
                 probe=row.names(r_probe),
                 codeword=row.names(r_codeword) ))
    
}
rep1_path <- 'path/to/rep1/out/directory'
rep1_path <- '/stornext/Bioinf/data/lab_phipson/data/xenium_human_breast_cancer/rep1/'

rep1 = get_xenium_data(rep1_path,
                mtx_name="cell_feature_matrix",
                trans_name="transcripts.csv.gz",
                cells_name="cells.csv.gz" )

graphclust_rep1=read.csv("graphclust_rep1.csv")
graphclust_rep1$anno =as.character(graphclust_rep1$Cluster)

t_cells =  c("CD4+_T_Cells","CD8+_T_Cells","T_Cell_&_Tumor_Hybrid",
             "Stromal_&_T_Cell_Hybrid")
dc_cells = c("LAMP3+_DCs","IRF7+_DCs")
macro_cells = c("Macrophages_1","Macrophages_2")
myo_cells = c("Myoepi_KRT15+", "Myoepi_ACTA2+")
tumor_cells = c("Invasive_Tumor", "Prolif_Invasive_Tumor")
graphclust_rep1[graphclust_rep1$Cluster %in% t_cells, "anno"] = "T_Cells"
graphclust_rep1[graphclust_rep1$Cluster %in% macro_cells,"anno"] = "Macrophages"
graphclust_rep1[graphclust_rep1$Cluster %in% c("DCIS_1", "DCIS_2"),"anno"] = "DCIS"
graphclust_rep1[graphclust_rep1$Cluster %in% dc_cells,"anno"] = "Dendritic"
graphclust_rep1[graphclust_rep1$Cluster %in% myo_cells,"anno"] = "Myoepithelial"
graphclust_rep1[graphclust_rep1$Cluster %in% tumor_cells,"anno"] = "Tumor"

graphclust_rep1$anno=factor(graphclust_rep1$anno,
                            levels=c("Tumor", "DCIS", "Stromal",
                                     "Macrophages","Myoepithelial",
                                     "T_Cells", "B_Cells",
                                     "Endothelial",
                                     "Dendritic", "Mast_Cells",
                                     "Perivascular-Like","Unlabeled"))


rep1_clusters = as.data.frame(cbind(as.character(graphclust_rep1$anno),
                                    paste("c",as.numeric(factor(graphclust_rep1$anno)),
                                          sep="")))

row.names(rep1_clusters) = graphclust_rep1$Barcode
colnames(rep1_clusters) = c("anno","cluster")

cells= rep1$cell_info
row.names(cells) = cells$cell_id
rp_names =  row.names(rep1_clusters)
rep1_clusters[rp_names,"x"] = cells[match(rp_names,row.names(cells)),
                                             "x_centroid"]
rep1_clusters[rp_names,"y"] = cells[match(rp_names,row.names(cells)),
                                             "y_centroid"]

rep1_clusters$anno=factor(rep1_clusters$anno,
                          levels=c("Tumor", "DCIS", "Stromal",
                                   "Macrophages","Myoepithelial",
                                   "T_Cells", "B_Cells",
                                   "Endothelial",
                                   "Dendritic", "Mast_Cells",
                                   "Perivascular-Like","Unlabeled"))
rep1_clusters$cluster=factor(rep1_clusters$cluster,
                             levels=paste("c", 1:12, sep=""))
table(rep1_clusters$cluster)
table(rep1_clusters$anno)
rep1_clusters$cells =row.names(rep1_clusters)
rep1_clusters$sample="rep1"

rep1_clusters = rep1_clusters[rep1_clusters$anno != "Unlabeled", ]


keep_cells=rep1_clusters[(rep1_clusters$x>=200 &
                              rep1_clusters$x<=500) &
                             (rep1_clusters$y>=200 &
                                  rep1_clusters$y<=1000),
                         "cells"]
subset_row_sums <- rowSums(rep1$cm[,keep_cells])

# Get the indices of the top 100 genes based on total expression
top_ids <- order(subset_row_sums, decreasing = TRUE)[1:50]

subset_genes = c("KRT7", "EPCAM", "FOXA1", # Tumor cells c1
                 "ERBB2", "SERPINA3", # DCIS cells c2
                 "LUM", "POSTN", "CCDC80",  # stromal cells  c3
                 "LYZ", "FGL2", "CD68",# Macrophages c4
                 "DST", "MYLK",  # Myoepithelial c5
                 "IL7R", "PTPRC", # T cell c6
                 "ITM2C", "MZB1", #B cells c7
                 "AQP1", "VWF", "PECAM1"  # Endothelial c8
                 )

#rm_genes = row.names(rep1$cm[,keep_cells])[rowSums(rep1$cm[,keep_cells]!=0) <=10]
rep1_sub_coord= rep1$trans_info[rep1$trans_info$cell_id %in% keep_cells &
                              rep1$trans_info$feature_name %in% subset_genes ,
                           c("x","y","feature_name","cell_id")]


rep1_clusters_sub =  rep1_clusters[rep1_clusters$cells %in% keep_cells,]
rep1_clusters_sub = rep1_clusters_sub[rep1_clusters_sub$cluster %in% paste("c", 1:8, sep=""),]
rep1_clusters_sub$cluster=factor(rep1_clusters_sub$cluster,
                                 levels=paste("c", 1:8, sep=""))




# create SpatialExperiment object
rep1_mol <- BumpyMatrix::splitAsBumpyMatrix(
    rep1_sub_coord[, c("x", "y")], 
    row = rep1_sub_coord$feature_name, col = rep1_sub_coord$cell_id )

rep1_sub<- SpatialExperiment(
    assays = list(molecules = rep1_mol),sample_id ="rep1" )

###################

# process the rep2 in the same way
rep2_path <- 'path/to/rep2/out/directory'
#rep2_path <- "/stornext/Bioinf/data/lab_phipson/data/xenium_human_breast_cancer/rep2/"
rep2 = get_xenium_data(rep2_path,
                mtx_name="cell_feature_matrix",
                trans_name="transcripts.csv.gz",
                cells_name="cells.csv.gz")

# with official clusters
graphclust_rep2=read.csv("graphclust_rep2.csv")
graphclust_rep2$anno = graphclust_rep2$Cluster
graphclust_rep2[graphclust_rep2$Cluster %in% t_cells, "anno"] = "T_Cells"
graphclust_rep2[graphclust_rep2$Cluster %in% macro_cells,"anno"] = "Macrophages"
graphclust_rep2[graphclust_rep2$Cluster %in% c("DCIS_1", "DCIS_2"),"anno"] = "DCIS"
graphclust_rep2[graphclust_rep2$Cluster %in% dc_cells,"anno"] = "Dendritic"
graphclust_rep2[graphclust_rep2$Cluster %in% myo_cells,"anno"] = "Myoepithelial"
graphclust_rep2[graphclust_rep2$Cluster %in% tumor_cells,"anno"] = "Tumor"

graphclust_rep2$anno=factor(graphclust_rep2$anno,
                            levels=c("Tumor", "DCIS", "Stromal",
                                     "Macrophages","Myoepithelial",
                                     "T_Cells", "B_Cells",
                                     "Endothelial",
                                     "Dendritic", "Mast_Cells",
                                     "Perivascular-Like","Unlabeled"))


rep2_clusters = as.data.frame(cbind(as.character(graphclust_rep2$anno),
                                    paste("c",as.numeric(factor(graphclust_rep2$anno)),
                                          sep="")))

row.names(rep2_clusters) = graphclust_rep2$Barcode
colnames(rep2_clusters) = c("anno","cluster")
rep2_clusters$cluster=factor(rep2_clusters$cluster)

cells= rep2$cell_info
row.names(cells) = cells$cell_id
rp_names =  row.names(rep2_clusters)
rep2_clusters[rp_names,"x"] = cells[match(rp_names,row.names(cells)),"x_centroid"]
rep2_clusters[rp_names,"y"] = cells[match(rp_names,row.names(cells)),"y_centroid"]

rep2_clusters$cluster=factor(rep2_clusters$cluster,
                             levels=paste("c",1:12, sep=""))
rep2_clusters$cells =row.names(rep2_clusters)
rep2_clusters$sample="rep2"
rep2_clusters = rep2_clusters[rep2_clusters$anno != "Unlabeled", ]


keep_cells2=rep2_clusters[(rep2_clusters$x>=200 & rep2_clusters$x<=500) &
                        (rep2_clusters$y>=200 & rep2_clusters$y<=1000),
                          "cells"]

rep2_sub_coord = rep2$trans_info[rep2$trans_info$cell_id %in% keep_cells2 &
                             rep2$trans_info$feature_name %in% subset_genes ,
                           c("x","y","feature_name","cell_id")]

rep2_clusters_sub =  rep2_clusters[rep2_clusters$cells %in% keep_cells2,]
rep2_clusters_sub = rep2_clusters_sub[rep2_clusters_sub$cluster %in% paste("c", 1:8, sep=""),]
rep2_clusters_sub$cluster=factor(rep2_clusters_sub$cluster,
                                 levels=paste("c", 1:8, sep=""))


# create SpatialExperiment object
rep2_mol <- BumpyMatrix::splitAsBumpyMatrix(
    rep2_sub_coord[, c("x", "y")], 
    row = rep2_sub_coord$feature_name, col = rep2_sub_coord$cell_id )

rep2_sub<- SpatialExperiment(
    assays = list(molecules = rep2_mol),sample_id ="rep2" )


########
# negative control detections
rep1_neg_coord = rep1$trans_info[rep1$trans_info$cell_id %in% keep_cells &
                                     (rep1$trans_info$feature_name %in% c(rep1$probe, rep1$codeword)) ,
                                 c("x","y","feature_name","cell_id")]

rep2_neg_coord = rep2$trans_info[rep2$trans_info$cell_id %in% keep_cells2 &
                                     (rep2$trans_info$feature_name %in% c(rep2$probe, rep2$codeword) ) ,
                                 c("x","y","feature_name","cell_id")]

shared_nc =  intersect(rep1_neg_coord$feature_name,rep2_neg_coord$feature_name)

rep1_added = cbind(x=NA, y=NA, 
                   feature_name=setdiff(rep2_neg_coord$feature_name,shared_nc),
                   cell_id = NA)
rep2_added = cbind(x=NA, y=NA, 
                   feature_name=setdiff(rep1_neg_coord$feature_name,shared_nc),
                   cell_id = NA)

rep1_neg_coord = rbind(rep1_neg_coord,rep1_added)
rep2_neg_coord = rbind(rep2_neg_coord,rep2_added)

rep1_neg_coord$category = "probe"
rep1_neg_coord[rep1_neg_coord$feature_name %in% c(rep1$codeword,rep2$codeword), "category"] = "codeword"

rep2_neg_coord$category = "probe"
rep2_neg_coord[rep2_neg_coord$feature_name %in% c(rep1$codeword,rep2$codeword), "category"] = "codeword"

rep1_neg_mol <- BumpyMatrix::splitAsBumpyMatrix(
    rep1_neg_coord[, c("x", "y","category")], 
    row = rep1_neg_coord$feature_name, col = rep1_neg_coord$cell_id )

rep1_neg<- SpatialExperiment(
    assays = list(molecules = rep1_neg_mol),sample_id ="rep1")

rep2_neg_mol <- BumpyMatrix::splitAsBumpyMatrix(
    rep2_neg_coord[, c("x", "y","category")], 
    row = rep2_neg_coord$feature_name, col = rep2_neg_coord$cell_id )

rep2_neg<- SpatialExperiment(
    assays = list(molecules = rep2_neg_mol),sample_id ="rep2")


rep1_clusters = rep1_clusters_sub
rep2_clusters = rep2_clusters_sub
# 
# saveRDS(rep1_sub, "rep1_sub.Rds")
# saveRDS(rep1_neg, "rep1_neg.Rds")
# write.csv(rep1_clusters, "rep1_clusters.csv")
# 
# saveRDS(rep2_sub, "rep2_sub.Rds")
# saveRDS(rep2_neg, "rep2_neg.Rds")
# write.csv(rep2_clusters, "rep2_clusters.csv")

# export rep1, rep1_clusters and rep1_neg
usethis::use_data(rep1_neg, overwrite = TRUE)
usethis::use_data(rep1_sub, overwrite = TRUE)
usethis::use_data(rep1_clusters, overwrite = TRUE)
# export rep2, rep2_clusters and rep2_neg
usethis::use_data(rep2_neg, overwrite = TRUE)
usethis::use_data(rep2_sub, overwrite = TRUE)
usethis::use_data(rep2_clusters, overwrite = TRUE)

