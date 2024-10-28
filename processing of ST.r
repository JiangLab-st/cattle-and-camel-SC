library(Seurat)
library(dplyr)
CreateS1000Object <- function(
  matrix_path,
  png_path,
  spot_radius = NULL,
  min.cells = 5,
  min.features = 100
  ){
  expr <- Seurat::Read10X(matrix_path, cell.column = 1)
  object <- Seurat::CreateSeuratObject(counts = expr,
                               assay = 'Spatial',
                               min.cells=min.cells,
                               min.features=min.features)
  #Image zoom rate
  cal_zoom_rate <- function(width, height){
    std_width = 1000
    std_height = std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
    if (std_width / std_height > width / height){
      scale = width / std_width
    }
    else{
      scale = height / std_height
    }
    return(scale)
  }
  #read png
  png <- png::readPNG(png_path)
  zoom_scale <-  cal_zoom_rate(dim(png)[2], dim(png)[1])
  #read barcode pos file
  ReadBarcodePos <- function(barcode_pos_path){
    barcode_pos <- read.table(gzfile(barcode_pos_path),header = F) %>%
      dplyr::rename(Barcode = V1 , pos_w = V2, pos_h = V3)
    return(barcode_pos)
  }
  #get barcode pos file path
  barcode_pos_path <- paste0(matrix_path,'/barcodes_pos.tsv.gz')
  barcode_pos <- ReadBarcodePos(barcode_pos_path = barcode_pos_path)
  barcode_pos <- barcode_pos %>% dplyr::filter(., Barcode %in% rownames(object@meta.data))
  #make spatial coord file for seurat S4 class
  coord <- data.frame(tissue = 1,
                      row = barcode_pos$pos_h,
                      col = barcode_pos$pos_w,
                      imagerow = barcode_pos$pos_h,
                      imagecol = barcode_pos$pos_w)
  rownames(coord) <- barcode_pos$Barcode
  #spot radius
  spot_radius_lib <- c(0.00063, 0.00179, 0.0027, 0.0039, 0.004, 0.0045, 0.005, NA, NA, NA, NA, NA, 0.0120)
  if(is.null(spot_radius)){
    spot_radius <- spot_radius_lib[as.numeric(gsub('L', '', strsplit(tail(strsplit(matrix_path, '/')[[1]],1), '_')[[1]][1]))]
  }else{
    spot_radius = spot_radius
  }
  #object
  sample1 <-  new(Class = "VisiumV1",
                  image = png,
                  scale.factors = Seurat::scalefactors(zoom_scale, 100, zoom_scale, zoom_scale),
                  coordinates = coord,
                  spot.radius = spot_radius,
                  assay = 'Spatial',
                  key = "sample1_")
  object@images <- list(sample1 = sample1)

  return(object)
}




source("/work/home/luo_funong/shitao/00.spaceX/03.result/kindey/09.script/CreateBmkObject.R")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)



object_L7 <- CreateS1000Object(matrix_path="/work/home/luo_funong/shitao/00.spaceX/03.result/kindey/05.AllheStat/BSTViewer_project/subdata/L7_heAuto",
          png_path="/work/home/luo_funong/shitao/00.spaceX/03.result/kindey/05.AllheStat/BSTViewer_project/he_roi_small.png", 
          min.cells =5, min.features = 100)

p2 <- SpatialFeaturePlot(object, features = "nCount_Spatial") + theme(legend.position = "right")

ggsave("1.pdf")

SpatialFeaturePlot(object_L7, features = c("DLK1","ACTA2","S100A4","TPM2","BTC","MYH11"))

object_L7_fillter_Rm <- object_L7[!grepl("^mt-", rownames(object_L7)), ]


object_L7_fillter_Rm_Nom <- SCTransform(object_L7_fillter_Rm, assay = "Spatial", verbose = TRUE, method = "poisson")
object_L7_fillter_Rm_Nom_PCA <- RunPCA(object_L7_fillter_Rm_Nom, assay = "SCT", verbose = FALSE)
object_L7_fillter_Rm_Nom_PCA <- FindNeighbors(object_L7_fillter_Rm_Nom_PCA, reduction = "pca", dims = 1:30)
object_L7_fillter_Rm_Nom_PCA <- FindClusters(object_L7_fillter_Rm_Nom_PCA, verbose = FALSE)
object_L7_fillter_Rm_Nom_PCA <- RunUMAP(object_L7_fillter_Rm_Nom_PCA,dims=1:20)
cols <- c( "#712143", "#b23469" , "#bd5870", "#e5726f"   , "#ec8d72", "#f1af81" , "#f6cc93", "#fae5a8", "#fdf7c2", "#f4f9c5", "#e3edb5", "#cae8b5" , "#a9dcba", "#83cfb7" , "#69aec5", "#6492c2" , "#7e74b4" , "#4f457f")
p1 <- DimPlot(object_L7_fillter_Rm_Nom_PCA, reduction = "umap", label = TRUE)
ggsave("sc_UMAP.pdf")

p2 <- SpatialDimPlot(object_L7_fillter_Rm_Nom_PCA, label = TRUE, label.size = 5)
p2 <- SpatialDimPlot(object_L7_fillter_Rm_Nom_PCA, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)

ggsave("spatial_UMAP.pdf")

de_markers <- FindAllMarkers(object_L7_fillter_Rm_Nom_PCA)
write.table(de_markers,"marker",sep="\t")
SpatialDimPlot(object_L7_fillter_Rm_Nom_PCA, cells.highlight = CellsByIdentities(object = object_L7_fillter_Rm_Nom_PCA, idents = c(9,10)), facet.highlight = TRUE, ncol = 3)
ggsave("spatial_UMAP_9_10.pdf")

de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
#brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
#    selection.method = "moransi")

#top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
#SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))


write.table(object_L4@images$sample1@coordinates, "split.table")
#L4_subset <- subset(object_L4, sample1_imagerow > 300 | sample1_imagecol < 300, invert = TRUE)
L4_subset <- subset(object_L4, sample1_imagerow > 700 | sample1_imagecol < 400, invert = TRUE)
p1 <- SpatialDimPlot(L4_subset, crop = TRUE, label = TRUE)
ggsave("L4_split.pdf")



Seurat_obj = readRDS('rds')
coor = data.frame('UMAP_1'= Seurat_obj@images$image@coordinates$imagerow , 'UMAP_2' = Seurat_obj@images$image@coordinates$imagecol)
rownames(coor) = colnames(Seurat_obj)
coor = as.matrix(coor)
Seurat_obj@reductions$umap@cell.embeddings = coor  
FeaturePlot(Seurat_obj, features = c('Pecam1','Cd38'), blend = TRUE)  
