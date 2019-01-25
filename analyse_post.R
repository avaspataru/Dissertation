library(Seurat)
library(dplyr)

pre.data <- read.delim("C:/Users/spata/Documents/Uni of warwick/cs3 - Project/post")

#remove genes column
genes <- pre.data$GENE
pre.data$GENE <- NULL

#add the genes column as row labels
row.names(pre.data)
row.names(pre.data) <- genes
.rowNamesDF(pre.data, make.names=FALSE) <- genes


# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
dataobj <- CreateSeuratObject(raw.data = pre.data, min.cells = 3, min.genes = 200,project = "NEWDATA2")

#STEP 2
mito.genes <- grep(pattern = "^MT-", x = rownames(x = dataobj@data), value = TRUE)
percent.mito <- Matrix::colSums(dataobj@raw.data[mito.genes, ])/Matrix::colSums(dataobj@raw.data)

#CellsMeta = data@data.info

dataobj <- AddMetaData(object = dataobj, metadata = percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))

dataobj <- FilterCells(object = dataobj, subset.names = c("percent.mito","nUMI","nGene"), 
                       low.thresholds = c(0.0,15000,0), high.thresholds = c(0.2,400000,9000))

dataobj <- NormalizeData(object = dataobj, normalization.method = "LogNormalize", 
                         scale.factor = 10000)

#VlnPlot(object = dataobj, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

#GenePlot(object = dataobj, gene1 = "nUMI", gene2 = "percent.mito")
#GenePlot(object = dataobj, gene1 = "nUMI", gene2 = "nGene")

dataobj <- FindVariableGenes(object = dataobj, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.2, x.high.cutoff = 3, y.cutoff = 0.5)
#length(x = dataobj@var.genes)

dataobj <- ScaleData(object = dataobj)

dataobj <- RunPCA(object = dataobj, pc.genes = dataobj@var.genes, 
                  do.print = TRUE, pcs.print = 1:100, genes.print = 5)

dataobj <- FindClusters(object = dataobj, reduction.type = "pca", dims.use = 1:10, 
                        resolution = 0.3, print.output = 0, save.SNN = TRUE, n.iter=5000, k.param=20)

PrintFindClustersParams(object = dataobj)
dataobj <- RunTSNE(object = dataobj, dims.use = 1:10, do.fast = TRUE,  check_duplicates = FALSE, n.iter=5000, perplexity=50)
TSNEPlot(object = dataobj)


# find markers for every cluster compared to all remaining cells, report
# only the positive ones
dataobj.markers <- FindAllMarkers(object = dataobj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
dataobj.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

top10 <- dataobj.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = dataobj, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

