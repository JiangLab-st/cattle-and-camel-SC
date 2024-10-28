library(Seurat)
library(ggplot2)
list <- read.table("./list")
list <- as.matrix(list)
list <- as.character(list)
files <- list.files(path="./", pattern='RDS$',full.names = TRUE)
files

objectlist <- list()
for(i in 1:length(files)){
   objectlist[[i]] <- readRDS(files[i])
 }
 
object.combined <- merge(x = objectlist[[1]], y = objectlist[2:3], add.cell.ids =list)
object.combined
object.combined <- NormalizeData(object.combined)
object.combined <- FindVariableFeatures(object.combined, selection.method = "vst", nfeatures = 2000)
object.combined <- ScaleData(object.combined, verbose = FALSE)
object.combined <- RunPCA(object.combined, npcs = 20, seed.use = 123,verbose = FALSE)
#cca
combined <- IntegrateLayers(object = object.combined,method = CCAIntegration,orig.reduction = "pca",new.reduction = "integrated.cca",verbose = FALSE)

combined <- FindNeighbors(combined, reduction = "integrated.cca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.1, cluster.name = "cca_clusters")
combined <- RunUMAP(combined, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")
DimPlot(combined,reduction = "umap.cca",group.by = "species",combine = FALSE,label.size = 2)
ggsave("12.pdf")

#harmony
combined <- IntegrateLayers(object = object.combined, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:20)
combined <- FindClusters(combined, resolution = 1, cluster.name = "harmony_clusters")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20, reduction.name = "umap")
DimPlot(combined,reduction = "umap",group.by = "sample",label.size = 2)
ggsave("12.pdf")

#rpca
combined <- IntegrateLayers(object = object.combined, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "integrated.rpca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.1, cluster.name = "rpca_clusters")
combined <- RunUMAP(combined, reduction = "integrated.rpca", dims = 1:20, reduction.name = "umap.rpca")
DimPlot(combined,reduction = "umap.rpca",group.by = "species",label.size = 2)
ggsave("12.pdf")
