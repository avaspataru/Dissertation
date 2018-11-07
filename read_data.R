library(Seurat)
library(dplyr)

rm(list = ls(all.names = TRUE))

pre.data <- read.delim("C:/Users/spata/Documents/Uni of warwick/cs3 - Project/pre.txt")

#remove genes column
genes <- pre.data$GENE
pre.data$GENE <- NULL

#add the genes column as row labels
row.names(pre.data)
row.names(pre.data) <- genes
.rowNamesDF(pre.data, make.names=FALSE) <- genes

#View(pre.data)

#post.data <- read.delim("C:/Users/spata/Documents/Uni of warwick/cs3 - Project/post")

#remove genes column
#genes1 <- post.data$GENE
#post.data$GENE <- NULL

#add the genes column as row labels
#row.names(post.data)
#row.names(post.data) <- genes1
#.rowNamesDF(post.data, make.names=FALSE) <- genes1
#View(post.data)

