library(Seurat)
library(dplyr)

rm(list = ls(all.names = TRUE))

new.data <- read.delim("C:/Users/spata/Documents/Uni of warwick/cs3 - Project/new dataset/file-DroNc4_N705_500.dge.txt")
#pre.data <- read.delim("C:/Users/spata/Documents/Uni of warwick/cs3 - Project/pre.txt")

#remove genes column
genes <- new.data$GENE
new.data$GENE <- NULL

#add the genes column as row labels
row.names(new.data)
row.names(new.data) <- genes
.rowNamesDF(new.data, make.names=FALSE) <- genes

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
#eliminate those expressed in < 1 cells
dataobj <- CreateSeuratObject(raw.data = new.data, min.cells = 1, min.genes = 200,project = "NEWDATA2")

#STEP 2
mito.genes <- grep(pattern = "^MT-", x = genes, value = TRUE)
percent.mito <- Matrix::colSums(dataobj@raw.data[mito.genes, ])/Matrix::colSums(dataobj@raw.data)

#CellsMeta = data@data.info

dataobj <- AddMetaData(object = dataobj, metadata = percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))

#dataobj <- FilterCells(object = dataobj, subset.names = c("percent.mito"), 
#                    low.thresholds = c(0.0), high.thresholds = c(0.05))

dataobj <- NormalizeData(object = dataobj, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

#VlnPlot(object = dataobj, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

#GenePlot(object = dataobj, gene1 = "nUMI", gene2 = "percent.mito")
#GenePlot(object = dataobj, gene1 = "nUMI", gene2 = "nGene")

dataobj <- FindVariableGenes(object = dataobj, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
#length(x = dataobj@var.genes)

dataobj <- ScaleData(object = dataobj)

dataobj <- RunPCA(object = dataobj, pc.genes = dataobj@var.genes, 
                  do.print = TRUE, pcs.print = 1:100, genes.print = 5)

dataobj <- FindClusters(object = dataobj, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = dataobj)
dataobj <- RunTSNE(object = dataobj, dims.use = 1:10, do.fast = TRUE,  check_duplicates = FALSE)
TSNEPlot(object = dataobj)