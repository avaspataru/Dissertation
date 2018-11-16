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

#genes are rows
#cells are columns

nouseless=0
for(i in names(new.data)){
  colum = new.data[,i]
  nogenes=0
  for(vari in colum){
    if(vari>0)
      nogenes=nogenes+1
  }
  if(nogenes<=200){
    new.data[[i]] <- NULL
    nouseless=nouseless+1
  }
}

print(nouseless)
print("taking out genes")

countg = 0
for(gene in genes){
  #print(new.data[gene,])
  reads=0
  for(entry in new.data[gene,]){
    #see if entry has more than 200 genes expressed
     if(entry>0)
          reads=reads+entry
  }
  if(reads == 0){
    countg=countg+1
  }
}

print(countg)

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
dataobj1 <- CreateSeuratObject(raw.data = new.data, min.cells = 1, min.genes = 200,project = "NEWDATA2")
