library(Seurat)
ss <- readRDS("liver.RDS")
gene1 <- list(c("ACADS", "MCEE", "MMAA", "MMUT", "PCCA", "PCCB", "PCK1", "PCK2", "PHYH"))
gene2 <- list(c("AACS", "ACAT1", "ACSS3", "BDH1", "BDH2", "HMGCL", "HMGCLL1", "HMGCS2", "SLC27A5"))
ss <- AddModuleScore( object = ss, features = gene1, name="SCFA_catabolic")
ss <- AddModuleScore( object = ss, features = gene2, name="ketone_body_biosynthetic")
Idents(ss) <- ss$tissue
table(Idents(ss))
#ss <- subset(x = ss, idents = c("Liver", "Kidney"), invert = TRUE)
VlnPlot(ss, features="SCFA_catabolic1",cols=c("#f59a9a", "#fe9533","#33b8b3", "#fde9a1","#747F7F"))
ggsave("SCFA_catabolic.pdf", width=12,height=12)
VlnPlot(ss, features="ketone_body_biosynthetic1",cols=c("#f59a9a", "#fe9533","#33b8b3", "#fde9a1","#747F7F"))
ggsave("ketone_body_biosynthetic.pdf", width=12,height=12)
#ggsave("SCFA_catabolic.pdf", width=12,height=12)
table <- as.data.frame(ss$SCFA_catabolic1)
write.table(table, "SCFA_catabolic", sep='\t')
table <- as.data.frame(ss$tissue)
write.table(table, "cell_tissue.table", sep='\t')
table <- as.data.frame(ss$ketone_body_biosynthetic1)
write.table(table, "ketone_body_biosynthetic", sep='\t')
#shell
paste cell_tissue.table ketone_body_biosynthetic|awk -v OFS='\t' '{print $2, $4}'|sed "1d" > ketone_body_biosynthetic.table
paste cell_tissue.table SCFA_catabolic|awk -v OFS='\t' '{print $2, $4}'|sed "1d" > SCFA_catabolic.table
vim

#R
#wilcox
df <- read.table("SCFA_catabolic.table",sep='\t', header=T)
with(data=df,
     pairwise.wilcox.test(x=score,g=tissue,p.adjust.method = "BH")
     )


Idents(ss) <- ss$tissue
table(Idents(ss))
AverageExp <- AverageExpression(ss)
expr <- AverageExp$RNA
write.table(expr, "expr.table", sep='\t')
saveRDS(ss, "Fibro_score_silei.RDS")
