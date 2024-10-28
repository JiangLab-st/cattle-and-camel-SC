library(Seurat)
ss <- readRDS("seurat_V5.RDS")
ss <- JoinLayers(ss,"RNA")
counts <- ss[["RNA"]]$counts 
ct <- as.matrix(counts) 
write.table(data.frame(rownames(ct),rownames(ct)),file = 'genes.tsv',
            quote = F,sep = '\t',
            col.names = F,row.names = F)
#barcodes.tsv
write.table(colnames(ct),file = 'barcodes.tsv',quote = F,
            col.names = F,row.names = F)
file="matrix.mtx"
sink(file)
cat("%%MatrixMarket matrix coordinate integer general\n")
cat("%\n")
cat(paste(nrow(ct),ncol(ct),sum(ct>0),"\n")) 
sink()
tmp=do.call(rbind,lapply(1:ncol(ct),function(i){
  return(data.frame(row=1:nrow(ct),
                    col=i,
                    exp=ct[,i]))
}) )
tmp=tmp[tmp$exp>0,]
head(tmp)
write.table(tmp,file = 'matrix.mtx',quote = F,
            col.names = F,row.names = F,append = T )


umap_XY <- ss@reductions$umap@cell.embeddings
head(umap_XY)
umap_XY <- as.data.frame(umap_XY)
write.table(umap_XY,"umap_XY",sep="\t")
write.table(ss@meta.data,"mata.table",sep="\t")
