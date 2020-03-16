library(garnett)
# data loading
mat <- Matrix::readMM(system.file("extdata", "exprs_sparse.mtx", package = "garnett"))
fdata <- read.table(system.file("extdata", "fdata.txt", package = "garnett"))
pdata <- read.table(system.file("extdata", "pdata.txt", package = "garnett"),
                    sep="\t")
row.names(mat) <- row.names(fdata)
colnames(mat) <- row.names(pdata)
str(mat)
length(mat[1,])
##create a new CDS object
pd <- new("AnnotatedDataFrame", data = pdata)
fd <- new("AnnotatedDataFrame", data = fdata)
pbmc_cds <- newCellDataSet(as(mat, "dgCMatrix"),
                           phenoData = pd,
                           featureData = fd)

#pbmc_cds_onlyMat <- newCellDataSet(as(mat, "dgCMatrix"))  #so the pd and fd are needed----> need start from the raw 10X data files

# generate size factors for normalization later
pbmc_cds <- estimateSizeFactors(pbmc_cds)

#------marker file----
#----marker QC
library(org.Hs.eg.db)
marker_file_path <- system.file("extdata", "pbmc_bad_markers.txt",
                                package = "garnett")
marker_check <- check_markers(pbmc_cds, marker_file_path,
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")
class(marker_check)# output data.frame
plot_markers(marker_check)
####################-------train the classifier---------##########
library(org.Hs.eg.db)
set.seed(260)

marker_file_path <- system.file("extdata", "pbmc_test.txt",
                                package = "garnett")
pbmc_classifier <- train_cell_classifier(cds = pbmc_cds,
                                         marker_file = marker_file_path,
                                         db=org.Hs.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 50,  #defaul= 500
                                         marker_file_gene_id_type = "SYMBOL")
head(pData(pbmc_cds))
class(pbmc_classifier)
###Viewing references
get_classifier_references(pbmc_classifier)

############-------------classify cell types#####################
pbmc_cds <- classify_cells(pbmc_cds, pbmc_classifier,
                           db = org.Hs.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "SYMBOL")

head(pData(pbmc_cds))

table(pData(pbmc_cds)$cell_type)
table(pData(pbmc_cds)$cluster_ext_type)
library(ggplot2)
qplot(tsne_1, tsne_2, color = cell_type, data = pData(pbmc_cds)) + theme_bw()
qplot(tsne_1, tsne_2, color = cluster_ext_type, data = pData(pbmc_cds)) + theme_bw()











#######################-------a second dataset PBMC2.7K-----------------#################
#### try PBMC2.7K
mat2<-Matrix::readMM("matrix.mtx")
fdata2 <- read.table("genes.tsv",sep = '\t')
class(fdata)
names(fdata)
row.names(fdata)

# modify fdata2
fdata2<-fdata2[,c(2,1)]
names(fdata2)<-c("gene_short_name", "geneID")
row.names(fdata2)<-fdata2$geneID
# modify pdata
pdata2 <- read.table("barcodes.tsv",sep="\t")
pdata2$Size_Factor<-1
row.names(pdata2)<-pdata2$V1
names(pdata2)[1]<-'cellID'

pd2 <- new("AnnotatedDataFrame", data = pdata2)
fd2 <- new("AnnotatedDataFrame", data = fdata2)
pbmc_cds2 <- newCellDataSet(as(mat2, "dgCMatrix"),
                            phenoData = pd2,
                            featureData = fd2)
pbmc_cds2 <- estimateSizeFactors(pbmc_cds2)

