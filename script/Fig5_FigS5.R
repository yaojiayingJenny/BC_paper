library(FlexDotPlot) #dotplot_cluster
library(tidyverse)
library(viridis)
library(ggalluvial)
library(ggsci)
library(gridExtra)
library(ComplexHeatmap)
library(Seurat)
require(dplyr)
require(Matrix)
require(magrittr)
library(scales)
library(ggplot2)
library(configr)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(harmony)
library(ggpubr) #ggboxplot




rds<-readRDS('rdsfiles/Epi_Minor.rds')
Epi_colors<-c('#D51F26','#208046','#2F2D66','#6E4B9E','#FBCB0A','#D24B27','#C06CAB','#D8A767','#89C75F','#90D5E4','#F37B7D','#9983BD','#8A9FD1')
prefix<-'Epi'
pdf(paste0(prefix,'_Minor_umap.pdf'),w=6.8,h=5)
p2 <- DimPlot(rds, reduction = "umap", label = T,group.by = "Minor",label.size = 3,pt.size = 0.3,raster=FALSE,cols=Epi_colors)+ggtitle('Epithelial cells')
print(p2)
dev.off()


rds<-readRDS('rdsfiles/Endo_Minor.rds')
Endo_colors<-c('#D51F26','#208046','#2F2D66','#6E4B9E','#FBCB0A','#D24B27','#C06CAB','#D8A767','#89C75F','#90D5E4','#F37B7D','#9983BD','#8A9FD1')
prefix<-'Endo'
pdf(paste0(prefix,'_Minor_umap.pdf'),w=6.8,h=5)
p2 <- DimPlot(rds, reduction = "umap", label = T,group.by = "Minor",label.size = 3,pt.size = 0.3,raster=FALSE,cols=Epi_colors)+ggtitle('Endothelial cells')
print(p2)
dev.off()

rds<-readRDS('rdsfiles/Fib_Minor.rds')
Fib_colors<-c('#D51F26','#208046','#2F2D66','#6E4B9E','#FBCB0A','#D24B27','#C06CAB','#D8A767','#89C75F','#90D5E4','#F37B7D','#9983BD','#8A9FD1')
prefix<-'Fib'
pdf(paste0(prefix,'_Minor_umap.pdf'),w=6.8,h=5)
p2 <- DimPlot(rds, reduction = "umap", label = T,group.by = "Minor",label.size = 3,pt.size = 0.3,raster=FALSE,cols=Epi_colors)+ggtitle('Fibroblast cells')
print(p2)
dev.off()

##---------------------------------
#FigS5 metabolism analysis
#https://github.com/wu-yc/scMetabolism
library(Seurat)
library(scMetabolism)
library(ggplot2)
library(rsvd)

prefix<-'BC'
rds<-readRDS('rdsfiles/BC_Minor.rds')
countexp.Seurat<-sc.metabolism.Seurat(obj = rds, method = "AUCell", imputation = F, ncores = 10, metabolism.type = "KEGG")


metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score #where metabolism.matrix
Barcodelist<-colnames(metabolism.matrix)
library(stringr)
Barcodelist1<-str_replace_all(Barcodelist,".1","-1")
colnames(metabolism.matrix)<-Barcodelist1
rds0<-CreateSeuratObject(counts = metabolism.matrix,min.cells=1,min.features=1)
rds0@meta.data<-rds@meta.data[rownames(rds0@meta.data),]
prefix<-'BC'
saveRDS(rds0,paste0(prefix,'_metabolism_last.rds'))

epi_colors<-c("#9ccc9c","#ff1100","#187bcd","#ff71b5")
tmp<-subset(rds0,celltype =='Epithelial')
# 1          Epi_IGFBP2 10753
# 2            Epi_CTGF  5376
# 3           Epi_UPK1A  2137
# 4   Epi_Proliferating   919
tmp$Minor<-factor(tmp$Minor,levels=c('Epi_IGFBP2','Epi_CTGF','Epi_UPK1A','Epi_Proliferating'))
Idents(tmp)<-'Minor'
table(tmp$Minor)
DefaultAssay(tmp) <- "RNA"
Idents(tmp)<-'Minor'
table(Idents(tmp))
cluster.averages <- AverageExpression(object = tmp, assays ='RNA',return.seurat = F)
write.csv(cluster.averages$RNA,file=paste(prefix,'Minor_averageMetabolism.csv',sep='_'),quote=F,row.names=T)
average.genes<-cluster.averages$RNA

tmp$Minor.g<-paste0(tmp$Minor,'.',tmp$Group)
Minor.glist<-c(paste0(c('Epi_IGFBP2','Epi_CTGF','Epi_UPK1A','Epi_Proliferating'),c('.Normal')),paste0(c('Epi_IGFBP2','Epi_CTGF','Epi_UPK1A','Epi_Proliferating'),c('.Tumor')))
tmp$Minor.g<-factor(tmp$Minor.g,levels=Minor.glist)
table(tmp$Minor.g)
Idents(tmp)<-'Minor.g'
table(Idents(tmp))
cluster.averages <- AverageExpression(object = tmp, assays ='RNA',return.seurat = F)
write.csv(cluster.averages$RNA,file=paste(prefix,'Minor.g_averageMetabolism.csv',sep='_'),quote=F,row.names=T)
average.genes<-cluster.averages$RNA

library(RColorBrewer)
colorlist<-rev(brewer.pal(n = 8, name = "RdBu"))
filename<-paste(prefix,'Minor.g_averageMetabolism_cluster-scale.pdf',sep='_')
pheatmap::pheatmap(average.genes,filename=filename,width=8,height=12,cluster_rows=T,cluster_cols=T,display_numbers = F,number_format = "%.0f",fontsize = 8,fontsize_col = 10,border_color=NA,angle_col = "45",color=color,scale='row')

saveRDS(tmp,file='Epi_metabolism_last.rds')







