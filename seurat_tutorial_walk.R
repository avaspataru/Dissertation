
seurat_tutorial <- function(given.data, project_name){
  
  # Initialize the Seurat object with the raw (non-normalized data).  Keep all
  # genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
  # least 200 detected genes
  dataobj <- CreateSeuratObject(raw.data = given.data, min.cells = 3, min.genes = 200,
                                project = "PREDATA")
  
  
  #STEP 2
  mito.genes <- grep(pattern = "^MT-", x = genes, value = TRUE)
  percent.mito <- Matrix::colSums(dataobj@raw.data[mito.genes, ])/Matrix::colSums(dataobj@raw.data)
  
  #CellsMeta = data@data.info
  
  data <- AddMetaData(object = data, metadata = percent.mito, col.name = "percent.mito")
  VlnPlot(object = data, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  
  
  # GenePlot is typically used to visualize gene-gene relationships, but can
  # be used for anything calculated by the object, i.e. columns in
  # object@meta.data, PC scores etc.  Since there is a rare subset of cells
  # with an outlier level of high mitochondrial percentage and also low UMI
  # content, we filter these as well
  par(mfrow = c(1, 2))
  GenePlot(object = data, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = data, gene1 = "nUMI", gene2 = "nGene")
}