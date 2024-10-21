# row seurat dat
# orthologous gene list
library(SingleCellExperiment)
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(harmony)
library(tidyr)
library(dplyr)
library(tibble)

setwd("xx")
dat <- c('Esophagus','Small_Intestine','Stomach','Salivary_Gland')
orth1 <- read.csv("orthologous gene list.csv", sep = ",", header = T)

RenameGenes <- function(obj,orth) {
  dat <- obj
  count <- dat@assays$RNA@counts
  count1 <- as.data.frame(count)
  count1$Cow.gene.name <- row.names(count1)
  orth2 <- orth1[,c("Gene.name","Cow.gene.name")]
  count1 <- left_join(count1,orth2,by='Cow.gene.name')
  count1 <- na.omit(count1)
  row.names(count1) <- count1$Gene.name
  count1 <- subset(count1, select = -c(Gene.name,Cow.gene.name))
  meta <- dat@meta.data
  obj <- CreateSeuratObject(count1,meta.data = meta)
  return(obj)
}

for (i in dat){
  cattle <- readRDS(paste0( i, '.RDS'))
  cattle1 <- RenameGenes(cattle,orth1)
  saveRDS(cattle1,paste0( i,'/cattle.RDS'))
}  
