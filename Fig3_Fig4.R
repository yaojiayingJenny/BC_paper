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



Mast<-readRDS('rdsfiles/Mast_Minor.rds')
Neu<-readRDS('rdsfiles/Neu_Minor.rds')


###--------------------------------------------------------------------
#Fig3a
celltype_colors<-c('#8BBDD3','#E06B4B','#9884AF','#A7592D')
prefix<-'Mast';rds<-Mast
pdf(paste0(prefix,'_Minor_umap.pdf'),w=6.8,h=5)
p2 <- DimPlot(rds, reduction = "umap", label = T,group.by = "Minor",label.size = 4,pt.size = 0.5,raster=FALSE,cols=celltype_colors)+ggtitle('Mast cells')
print(p2)
dev.off()



###--------------------------------------------------------------------
#Fig3f

celltype_colors<-c('#3C80B0','#379548','#ED8F39')
prefix<-'Neu';rds<-Neu
pdf(paste0(prefix,'_Minor_umap.pdf'),w=6.8,h=5)
p2 <- DimPlot(rds, reduction = "umap", label = T,group.by = "Minor",label.size = 4,pt.size = 0.5,raster=FALSE,cols=celltype_colors)+ggtitle('Neutrophil cells')
print(p2)
dev.off()


###--------------------------------------------------------------------
#Fig4a
rds<-readRDS('rdsfiles/TNK_Minor.rds')
celltype_colors<-c("#DD626F","#388BC0","#C4E7BF","#66C2A5","#F46D43","#FDAE61","#E6298A","#B15928","#6D419C","#CAB3D6","#426E66","#FED976","#225EA8","#F768A1","#41AB5D")
prefix<-'TNK'
pdf(paste0(prefix,'_Minor_umap.pdf'),w=6.8,h=5)
p2 <- DimPlot(rds, reduction = "umap", label = T,group.by = "Minor",label.size = 3,pt.size = 0,raster=FALSE,cols=celltype_colors)+ggtitle('T&NK cells')
print(p2)
dev.off()

###--------------------------------------------------------------------
#FigS4c
tmp<-rds;g1<-'NK_FGFBP2';g2<-'NK_XCL2'
type<-'NK';vs<-'_NK_FGFBP2_vs_NK_XCL2'
Idents(tmp)<-'Minor'
genes <- FindMarkers(tmp, ident.1 = g1, ident.2 =g2 ,logfc.threshold =0,test.use ='wilcox' ,min.pct = 0.1,min.cells.group=1)
genes$gene<-rownames(genes)
genes<-arrange(genes,desc(avg_log2FC))
diffgene<-subset(genes,p_val_adj<0.05 & ((avg_log2FC>=0.25 & pct.1>=0.1) |(avg_log2FC<=(-0.25) & pct.2>=0.1)))
diffgene$group <- ifelse(diffgene$avg_log2FC>0,'up','down')
'%!in%' <- function(x,y)!('%in%'(x,y))
genes1<-subset(genes,gene %!in%  diffgene$gene)
genes1$group<-'nonDE'
genes<-rbind(genes1,diffgene)
genes<-arrange(genes,desc(avg_log2FC))
write.csv(diffgene,paste(type,vs,"diffgene_FDR0.05.csv",sep=""),quote=F)
write.csv(genes,paste(type,vs,"diffgene_all.csv",sep=""),quote=F)
require("ggrepel") #geom_text_repel
#data$group <- ifelse(data$order<=20,'yes','no')
de<-c(head(diffgene,10)$gene,tail(diffgene,10)$gene)
data<-genes
main<-paste0(type,',',g1,' vs ',g2)
p0<-theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
legend.title = element_blank(),legend.text = element_text(face="plain", color="black",size=10),
axis.text.x = element_text(color="black",size=10,angle=0,hjust=0),
axis.text.y = element_text(color="black",size=10),
axis.title.x = element_text(face="plain", color="black",size=12),
axis.title.y = element_text(face="plain", color="black",size=12))+theme(legend.position='none')
p2<-ggplot(data,aes(x=avg_log2FC,y=-log10(p_val_adj),color=group))+ geom_point(shape=20)+theme_bw()+p0+geom_text_repel(data=data[de,], aes(label=gene),size=2,direction="both",min.segment.length = 0.05,segment.alpha=0.6,label.padding = 0.4,max.overlaps =30)+ggtitle(main)+scale_color_manual(values=c('#056FBB','#808080','#CE2314'),name="DE genes")+theme(legend.position= 'right')+geom_vline(xintercept=0,linetype='dashed')
pdf(paste0(type,vs,"diffgene_volcanic.pdf"),w=6.5,h=5)
print(p2)
dev.off()



