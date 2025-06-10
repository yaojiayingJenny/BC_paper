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



rds<-readRDS('rdsfiles/Mye_Minor.rds')


###--------------------------------------------------------------------
#1.Fig2a umap
prefix<-'Mye'
Mye_colors<-c('#D51F26','#208046','#2F2D66','#6E4B9E','#FBCB0A','#D24B27','#C06CAB','#D8A767','#89C75F','#90D5E4','#F37B7D','#9983BD','#8A9FD1')
pdf(paste0(prefix,'_Minor_umap.pdf'),w=6.8,h=5)
p2 <- DimPlot(rds, reduction = "umap", label = T,group.by = "Minor",label.size = 3,pt.size = 0.3,raster=FALSE,cols=Mye_colors)+ggtitle('Myeloid cells')
print(p2)
dev.off()

###--------------------------------------------------------------------
#2.FigSa rename genes dotplot
rename<-read.csv('Fig2/Mye_rename.xls',header=T,sep='\t')

# Minor	Markergenes	Major
# Mono1_CD14	S100A4,CD14,S100A8,S100A9,NLRP3	Mono
# Mono2_CD14	S100A4,CD14,S100A8,S100A9,S100A12	Mono
# Mono_CD16	S100A4,FCGR3A,CDKN1C,LILRB2	Mono
# Mac_TREM2	CD68,C1QA,C1QC,APOE,GPNMB,TREM2	Mac
# Mac_CCL3	CD68,C1QA,C1QC,CCL2,CCL3,CCL4	Mac
# Mac_IL1B	CD68,C1QA,C1QC,IL1B,CCL20,CXCL2,CXCL3	Mac
# Mac_CD163	CD68,C1QA,C1QC,MAF,F13A1,CD163,LILRB5,MARCO	Mac
# cDC_CD1C	FCER1A,CD1C	DC
# cDC_CD207	CD207,CD1A,CD1E	DC
# cDC_CLEC9A	CLEC9A,WDFY4,IRF8	DC
# DC_LAMP3	LAMP3,CCR7	DC
# Mye_IRF7	IRF7	DC
# Mye_Proliferating	MKI67,TUBB,TOP2A	MKI67


prefix<-'Mye'
tmp<-rds
genes<-unique(unlist(strsplit(paste(rename$Markergenes,collapse=","),split = ",",fixed=T)))
setdiff(genes,rownames(tmp))
tmp$Minor<-factor(tmp$Minor,levels=rev(map$V2))
Idents(tmp)<-'Minor'
tmp <- ScaleData(tmp, features = genes, verbose = FALSE)
p0=DotPlot(tmp,features =genes,group.by = "Minor")#,'FTH1','FTL'
p3<-p0+theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+xlab('')+ylab('')+
  guides(color = guide_colorbar(title = 'Scaled expression',order = 1),size = guide_legend("Percent expressed"),order = 0)+
  scale_colour_gradientn(colours = c("dodgerblue1","#44c1f0", "lightgoldenrod","#e20612",'#cc340c'))+
  theme(legend.position = "top")+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.3, linetype="solid"))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
pdf(paste0(prefix,'_rename_genes_dotplot_new.pdf'),w=10,h=6,onefile=F)
print(p3)
dev.off()


###--------------------------------------------------------------------
#3.Fig2b  #signature score vlnplot

# Antigen presenting/co-stimulatory capacity 
# M1 signature
# M2 signature

genes1<-read.csv('Fig2/Antigen_presenting_genelist.xls',header=T)[,1]
genes2<-read.csv('Fig2/M1_genelist.xls',header=T)[,1]
genes3<-read.csv('Fig2/M2_genelist.xls',header=T)[,1]

genes1
setdiff(genes1,rownames(Mye))
genes2
setdiff(genes2,rownames(Mye))
genes3
setdiff(genes3,rownames(Mye))
#modulescore
rds <- AddModuleScore(object = rds,features = list(genes1),name = 'Antigen')
rds <- AddModuleScore(object = rds,features = list(genes2),name = 'M1')
rds <- AddModuleScore(object = rds,features = list(genes3),name = 'M2')

Idents(rds)<-'Minor'
prefix<-'Antigen'
p<-VlnPlot(rds, features = 'Antigen1',pt.size = -1,slot="data",assay="RNA",stack = F,flip =TRUE,cols=Mye_colors,group.by='Minor')+theme(legend.position='none',axis.text.y = element_text(color="black",size=8,angle=0),axis.text.x = element_text(color="black",size=8,angle=30),axis.title.x=element_text(size = 8),axis.title.y=element_text(size = 10))+xlab('')+ylab('ModuleScore')+ggtitle('Antigen presenting/co-stimulatory capacity')+stat_summary (fun.data = function (x) data.frame (y=max(rds$Antigen1)-0.1, label = paste (round(mean (x),2))), geom="text")
pdf(paste0('Mye_celltype_',prefix,'_modulescore_VlnPlot.pdf'),w=8,h=2.5)
print(p)
dev.off()


prefix<-'M1'
p<-VlnPlot(rds, features = 'M11',pt.size = -1,slot="data",assay="RNA",stack = F,flip =TRUE,cols=Mye_colors,group.by='Minor')+theme(legend.position='none',axis.text.y = element_text(color="black",size=8,angle=0),axis.text.x = element_text(color="black",size=8,angle=30),axis.title.x=element_text(size = 8),axis.title.y=element_text(size = 10))+xlab('')+ylab('ModuleScore')+ggtitle('M1 signature')+stat_summary (fun.data = function (x) data.frame (y=max(rds$M11)-0.1, label = paste (round(mean (x),2))), geom="text")
pdf(paste0('Mye_celltype_',prefix,'_modulescore_VlnPlot.pdf'),w=8,h=2.5)
print(p)
dev.off()


prefix<-'M2'
p<-VlnPlot(rds, features = 'M21',pt.size = -1,slot="data",assay="RNA",stack = F,flip =TRUE,cols=Mye_colors,group.by='Minor')+theme(legend.position='none',axis.text.y = element_text(color="black",size=8,angle=0),axis.text.x = element_text(color="black",size=8,angle=30),axis.title.x=element_text(size = 8),axis.title.y=element_text(size = 10))+xlab('')+ylab('ModuleScore')+ggtitle('M2 signature')+stat_summary (fun.data = function (x) data.frame (y=max(rds$M21)-0.1, label = paste (round(mean (x),2))), geom="text")
pdf(paste0('Mye_celltype_',prefix,'_modulescore_VlnPlot.pdf'),w=8,h=2.5)
print(p)
dev.off()



