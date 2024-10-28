#Rscript  Cyto_Lactation.R -I ./mtx/  -M ./meta.table -XY ./UMAP_XY  -G TP63
library(Seurat)
library(CytoTRACE)
library(dplyr)
library(ggplot2)
Args <- commandArgs()

#Get the parameters
parser = argparse::ArgumentParser(description="Cyte")
parser$add_argument('-I','--input', help='input mtx dir with 3 files')
parser$add_argument('-M','--meta', help='input original metatable')
parser$add_argument('-X','--XY', help='original UMAP_XY')
parser$add_argument('-G','--gene',help='genename')
#parser$add_argument('-PC','--pc',help='pc usage')

args = parser$parse_args()

EC.data <- Read10X(data.dir = args$input, gene.column = 1)
meta <- read.table(args$meta,header=T,sep='\t')
#ifnb.data <- readRDS(args$input)
#table(ifnb.data$orig.ident)
#table(ifnb.data$celltype)
ifnb.data <- CreateSeuratObject(EC.data,meta.data = meta)
object.combined <- NormalizeData(ifnb.data)
object.combined <- FindVariableFeatures(object.combined, selection.method = "vst", nfeatures = 2000)
object.combined <- ScaleData(object.combined, verbose = FALSE)
object.combined <- RunPCA(object.combined, npcs = 20, seed.use = 123,verbose = FALSE)
ifnb.data <- RunUMAP(object.combined, reduction = "pca", dims = 1:20, seed.use = 123)
umap_XY <- read.table(args$XY,header=T,sep='\t')
colnames(umap_XY) <- c('UMAP_1','UMAP_2')
loca <- as.matrix(umap_XY)
a <- match(rownames(ifnb.data@meta.data),rownames(loca))
loca <- loca[a,]
ifnb.data@reductions[["umap"]]@cell.embeddings[,1] <- loca[,1]
ifnb.data@reductions[["umap"]]@cell.embeddings[,2] <- loca[,2]
saveRDS(ifnb.data, "new_V4_ori.UMAPmeta.RDS")
