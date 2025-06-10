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


#Fig5a
rds<-readRDS('rdsfiles/Fib_Minor.rds')
Fib_colors<-c('#D51F26','#208046','#2F2D66','#6E4B9E','#FBCB0A','#D24B27','#C06CAB','#D8A767','#89C75F','#90D5E4','#F37B7D','#9983BD','#8A9FD1')
prefix<-'Fib'
pdf(paste0(prefix,'_Minor_umap.pdf'),w=6.8,h=5)
p2 <- DimPlot(rds, reduction = "umap", label = T,group.by = "Minor",label.size = 3,pt.size = 0.3,raster=FALSE,cols=Epi_colors)+ggtitle('Fibroblast cells')
print(p2)
dev.off()


#Fig5e
rds<-readRDS('rdsfiles/Epi_Minor.rds')
Epi_colors<-c('#D51F26','#208046','#2F2D66','#6E4B9E','#FBCB0A','#D24B27','#C06CAB','#D8A767','#89C75F','#90D5E4','#F37B7D','#9983BD','#8A9FD1')
prefix<-'Epi'
pdf(paste0(prefix,'_Minor_umap.pdf'),w=6.8,h=5)
p2 <- DimPlot(rds, reduction = "umap", label = T,group.by = "Minor",label.size = 3,pt.size = 0.3,raster=FALSE,cols=Epi_colors)+ggtitle('Epithelial cells')
print(p2)
dev.off()

#FigS6c
rds<-readRDS('rdsfiles/Endo_Minor.rds')
Endo_colors<-c('#D51F26','#208046','#2F2D66','#6E4B9E','#FBCB0A','#D24B27','#C06CAB','#D8A767','#89C75F','#90D5E4','#F37B7D','#9983BD','#8A9FD1')
prefix<-'Endo'
pdf(paste0(prefix,'_Minor_umap.pdf'),w=6.8,h=5)
p2 <- DimPlot(rds, reduction = "umap", label = T,group.by = "Minor",label.size = 3,pt.size = 0.3,raster=FALSE,cols=Epi_colors)+ggtitle('Endothelial cells')
print(p2)
dev.off()

##---------------------------------
#FigS7f metabolism analysis
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



#FigS7f绘图
#ptest
ge1<-'Glycolysis / Gluconeogenesis';ge2<-'Glycolysis-Gluconeogenesis'
ge1<-'Citrate cycle (TCA cycle)';ge2<-'Citrate cycle (TCA cycle)'
dataset<-FetchData(tmp,vars=ge1)
meta<-tmp@meta.data[rownames(dataset),c('Minor','Group')]
dataset<-cbind(dataset,meta)
#
library(ggsignif)	
ptest1<-stat_compare_means(comparisons =list(c('Normal','Tumor')),method='wilcox.test',paired=FALSE,size=5,vjust = 0.5,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),symbols = c("***", "**", "*", "ns")),hide.ns = TRUE)

#label = "p.format"
#label = "p.signif"
#ptest1<-stat_compare_means(label = "p.signif",comparisons =list(c('AC Bln','AC Ag'),c('AA Bln','AA Ag')),method='wilcox.test',paired=FALSE,size=5,vjust = 0.5,hide.ns = TRUE)

p0<-theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
		legend.title = element_blank(),legend.text = element_text(face="plain", color="black",size=15),
		axis.text.x = element_text(color="black",size=10,angle=30,hjust=1,vjust=1),
		axis.text.y = element_text(color="black",size=10),
		axis.title.x = element_text(face="plain", color="black",size=15),
		axis.title.y = element_text(face="plain", color="black",size=10))+theme(legend.position='none')
dataset$Minor<-factor(dataset$Minor,levels=c('Epi_IGFBP2','Epi_CTGF','Epi_UPK1A','Epi_Proliferating'))	
#genelist<-'CCL17,CXCL16,TREM2,GPNMB'
dataset$value<-dataset[,1]
pair_groups4<-c('#69B38D','#339191','#E19056','#E6612C')
pair_groups4<-c('#6BB7CA','#ED9F9B')
p<-ggplot(dataset, aes(x = Group, y = value,fill=Group)) +geom_violin()+facet_wrap(~Minor,ncol=4)+ptest1+theme_bw()+p0+theme(legend.position='none')+xlab('')+scale_fill_manual(values=pair_groups4)+ylim(0,max(dataset$value)+0.02)+ylab(colnames(dataset)[1])+stat_summary (fun.data = function (x) data.frame (y=max(dataset$value), label = paste (round(mean (x),4))), geom="text",size=2)
pdf(paste0(prefix,'_',ge2,'_celltype_group_VlnPlot_wilcox-test.pdf'),w=8,h=2.5)
print(p)
dev.off()




