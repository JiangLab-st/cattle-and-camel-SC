#Get preliminary alignment results from eggnog specie.csv
#Get annotation information from gtf file species.txt
library(dplyr)
library(stringr)
#Protein stable ID version
cattle1 <- read.csv('xxx/catttle.CSV',sep=',',header=T)
cattle2 <- read.csv('xxx/cattle.txt',sep = '\t',header = T)
cattle <- left_join(cattle1,cattle2,by='Protein.stable.ID.version')
cattle <- cattle[,c(1,9,22,23,24,25)]
cattle1 <- cattle[cattle$Gene.name=="",]
cattle1$Gene.name <- cattle1$Gene.stable.ID
cattle2 <- cattle[!cattle$Gene.name=="",]
cattle <- rbind(cattle1,cattle2)
cattle$gene <- paste0(cattle$Gene.stable.ID,'_',cattle$Gene.start..bp.,'_',cattle$Gene.end..bp.)
cattle1 <- cattle[!duplicated(cattle$gene), ]
cattle1$Preferred_name <- toupper(cattle1$Preferred_name)
dup_rows <- duplicated(cattle1$Preferred_name) | duplicated(cattle1$Preferred_name, fromLast = TRUE)
cattle1$perfectname_X <- ifelse(dup_rows, 1, 0)
dup_rows <- duplicated(cattle1$Gene.name) | duplicated(cattle1$Gene.name, fromLast = TRUE)
cattle1$Gene.name_X <- ifelse(dup_rows, 1, 0)
cattle1 <- cattle1[cattle1$Gene.name_X==0,]
cattle1 <- cattle1[cattle1$perfectname_X==0,]
write.csv(cattle1,'cattle1.csv')

#Integrate orthologous genes from multiple species
cattle <- read.csv('xxx/catttle1.csv',sep=',',header=T)  #The result of the previous step
camel <- read.csv('xxx/camel1.csv',sep=',',header=T)
dat <- left_join(cattle,camel,by='Preferred_name') #Multiple species are added sequentially
old_column_names <- c("Preferred_name", "Gene.name_cow", "Gene.name_camel")
new_column_names <- c("Gene.name", "Cow.gene.name", "Camel.gene.name")
for (i in 1:length(new_column_names)) {
  colnames(dat)[colnames(er) == old_column_names[i]] <- new_column_names[i]
}
write.csv(dat,'xxx/cattle-camel-orthologous.CSV')
