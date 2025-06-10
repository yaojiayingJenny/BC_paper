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

rds<-readRDS('rdsfiles/BC_celltype.rds')

###--------------------------------------------------------------------
#1.Fig1c umap
celltype_colors<-read.csv('Fig1/celltype_colors.xls',check.names=F,header=F)$V1
# E3879E
# 6DA34D
# 9D4EDD
# 6AB8EE
# F433AB
# F44336
# 7180B9
# B68D40
# F8961E
b<-DimPlot(rds, group.by = "celltype", label = TRUE, repel = TRUE,pt.size = 0,raster=FALSE,cols=celltype_colors) + NoLegend() + xlim(-10, 15) + ylim(-16, 11)+ labs(x = "UMAP1", y = "UMAP2") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+ggtitle('')
pdf('Fig1c.umap_all_celltype.pdf',w=5,h=5);print(b);dev.off()


###--------------------------------------------------------------------
#2.Fig1f cell proporatio test
prefix<-'BC'
tmp0<-rds
a<-table(tmp0$Sample,tmp0$celltype)
write.csv(a,file=paste0(prefix,'_celltype_sample_cellnum.csv'),quote=F)
b<-data.frame(table(tmp0@meta.data[,c('Sample','celltype')]),margin=1)
per<-round(prop.table(table(tmp0@meta.data[,c('Sample','celltype')]), margin=1),3)
per <-data.frame(per)
per$cellnum<-b$Freq
dataset<-read.table('group_info.xls',header=T,sep='\t')
#group_info.xls
# PatientID	SampleID	Group
# BC1	BC1P	Normal
# BC2	BC2P	Normal
# BC3	BC3P	Normal
# BC4	BC4P	Normal
# BC5	BC5P	Normal
# BC6	BC6P	Normal
# BC7	BC7P	Normal
# BC8	BC8P	Normal
# BC9	BC9P	Normal
# BC10	BC10P	Normal
# BC11	BC11P	Normal
# BC1	BC1T	Tumor
# BC2	BC2T	Tumor
# BC3	BC3T	Tumor
# BC4	BC4T	Tumor
# BC5	BC5T	Tumor
# BC6	BC6T	Tumor
# BC7	BC7T	Tumor
# BC8	BC8T	Tumor
# BC9	BC9T	Tumor
# BC10	BC10T	Tumor
# BC11	BC11T	Tumor

per$Group<- plyr::mapvalues(x =per$Sample ,from =dataset[,2],to = dataset[,3])
print("每个样本在所有细胞簇中的占比情况如下：")
head(per)
write.table(per,file=paste0(prefix,'_celltype_sample_group_ratio.xls'),quote=F,row.names=F,sep='\t')
Group_colors<-c("#4DBBD5FF","#E64B35FF")

rename<-read.csv('Fig1/celltype_rename.xls',header=T,sep='\t')
#Fig1/celltype_rename.xls
# celltype	Markergenes
# T&NK	PTPRC,CD3D,GNLY
# B	PTPRC,CD79A,MS4A1
# Plasma	PTPRC,CD79A,MZB1
# MonoMøDC	PTPRC,CSF1R,CD68,CD14,FCGR3A,CD1C
# Neutrophil	PTPRC,S100A8,S100A9,CSF3R
# Mast	PTPRC,MS4A2,TPSAB1,TPSB2,CPA3
# Epithelial	EPCAM,KRT19
# Endothelial	PECAM1,VWF
# Fibroblast	COL1A1,FN1

#FigS1c
per$celltype<-factor(per$celltype,levels=c(rename$celltype))
library(ggpubr) #ggboxplot
p0<-theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
        legend.title = element_blank(),legend.text = element_text(face="plain", color="black",size=15),
        axis.text.x = element_text(color="black",size=10,angle=30,hjust=0.5),
        axis.text.y = element_text(color="black",size=10),
        axis.title.x = element_text(face="plain", color="black",size=12),
        axis.title.y = element_text(face="plain", color="black",size=12))+theme(legend.position='none')
ptest<-stat_compare_means(comparisons =list(c('Normal','Tumor')),method='t.test',paired=TRUE,size=3,vjust = 1,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),symbols = c("***", "**", "*", "ns")))
p<-ggboxplot(per, 'Group', "Freq",
          color = 'Group', add = "jitter",size=0.3,
          xlab = " ", ylab = 'Cell proportion',palette = Group_colors) + rotate_x_text()+p0+facet_wrap(~celltype,ncol=5,scale='free')+ptest
pdf(paste(prefix,"celltype_group_cellratio_boxplot_paired-t.test.pdf",sep='_'),w=10,h=4)
print(p)
dev.off()
write.table(per,file=paste0(prefix,'_celltype_sample_group_ratio.xls'),quote=F,row.names=F,sep='\t')


#perform in each cells
celllist<-rename$celltype
for (type in celllist){
da<-subset(per,celltype==type)
p<-ggpaired(da, x = "Group", y = "Freq",
         color = "Group", line.color = "#808080", line.size = 0.4,linetype ='dashed',
         palette = Group_colors,xlab = " ", ylab = 'Cell proportion')+ rotate_x_text()+p0+ggtitle(type)+stat_compare_means(comparisons =list(c('Normal','Tumor')),method='t.test',paired=TRUE,size=3,hjust = 0.1,vjust = 0.1)
pdf(paste0(type,'_group_cellratio_boxplot_paired-t.test.pdf'),w=3,h=4)
print(p)
dev.off()
}

###--------------------------------------------------------------------
#Fig1d rename genes dotplot 
genes<-unique(unlist(strsplit(paste(rename$Markergenes,collapse=","),split = ",",fixed=T)))
setdiff(genes,rownames(tmp))
tmp$celltype<-factor(tmp$celltype,levels=rev(unique(rename$celltype)))
Idents(tmp)<-'celltype'
table(tmp$celltype)
tmp <- ScaleData(tmp, features = genes, verbose = FALSE)
p0=DotPlot(tmp,features =genes,group.by = "celltype")#,'FTH1','FTL'
p3<-p0+theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+xlab('')+ylab('')+
  guides(color = guide_colorbar(title = 'Scaled expression',order = 1),size = guide_legend("Percent expressed"),order = 0)+
  scale_colour_gradientn(colours = c("dodgerblue1","#44c1f0", "lightgoldenrod","#e20612",'#cc340c'))+
  theme(legend.position = "top")+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.3, linetype="solid"))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
pdf(paste0(prefix,'_celltype_rename_genes_dotplot.pdf'),w=8,h=5,onefile=F)
print(p3)
dev.off()


###--------------------------------------------------------------------
#Fig1e celltype ratio in each sample
prefix<-'BC'
a<-table(rds$celltype,rds$Sample)
write.csv(a,file=paste0(prefix,'_celltype_sample_cellnum.csv'),quote=F)
ratio<-round(table(rds$Sample,rds$celltype)/as.vector(table(rds$Sample)),4)
write.csv(ratio,file=paste0(prefix,'_celltype_sample_cellratio.csv'),quote=F)
ratio<-as.matrix(ratio)
a<-data.frame(CellType=rds@meta.data[,'celltype'],Samples=rds@meta.data[,'Sample'])
a$CellType<-factor(a$CellType,levels=sort(unique(a$CellType)))

p0<-theme(legend.position='right',panel.grid=element_blank(), legend.background = element_rect(colour = NA),
		legend.title = element_blank(),legend.text = element_text(face="plain", color="black",size=8),
		axis.text.x = element_text(color="black",size=10,angle=90, hjust = 1),
		axis.text.y = element_text(color="black",size=10),
		axis.title.x = element_text(face="plain", color="black",size=12),
		axis.title.y = element_text(face="plain", color="black",size=12))
pdf(paste(prefix,"celltype_sample_cellratio_barplot.pdf",sep='_'),w=10,h=5)
p<-ggplot(a, aes(Samples)) + geom_bar(aes(fill=CellType), position='fill',width=0.6)+labs(x=" ", y = "Cell type distribution",fill= "CellType")+theme_bw()+p0+scale_fill_manual(values=celltype_colors)
print(p)
dev.off()



##---------------------------------------
#python code
import scanpy
import scanpy as sc
import pandas as pd
import os 

#sc._settings.ScanpyConfig._vector_friendly = True
scanpy._settings.ScanpyConfig(_vector_friendly=True)
adata = sc.read_h5ad('Fig1/BC_Minor_seurat.rds.h5ad')


fig = ax = sc.pl.embedding_density(
    adata, basis='umap', groupby='Group', ncols=2, 
    color_map='plasma', return_fig=False,  #plasma color
    save="FigS1a.Group_density_umap.pdf")


fig = ax = sc.pl.embedding_density(
    adata, basis='umap', groupby='Patient', ncols=6, 
    color_map='plasma', return_fig=False, 
    save="FigS1b.Patient_density_umap.pdf")


#########
