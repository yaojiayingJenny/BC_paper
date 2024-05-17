library(FlexDotPlot) #dotplot_cluster
library(Seurat)
library(dplyr)
library(tidyverse)
library(viridis)
library(ggalluvial)
library(ggsci)
library(ggplot2)
library(gridExtra)
library(ComplexHeatmap)
require(Matrix)
require(magrittr)
library(scales)
library(configr)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(harmony)
library(SeuratData)
library(ggpubr) #ggboxplot

#1.QC 
```
We used DoubletDetection to detect double cells in each sample,respectively, and removed them in the following analysis.
min.cells = 3
Gene number: 200-10000(min_nfeature_rna = 200,max_nfeature_rna = 10000)
mt<=15%,HB<=5
```




#2.Integrate and cluster
prefix<-'BC';nfeatures<-2000;npc<-30;dims<-20
tmp0 <- NormalizeData(tmp0) %>% FindVariableFeatures(nfeatures =nfeatures) %>% ScaleData() %>% RunPCA(npcs = as.numeric(npc),verbose = FALSE)
tmp0 <- RunHarmony(tmp0, group.by.vars = "stim")
tmp0 <- RunUMAP(tmp0, reduction = "harmony", dims = 1:dims)
tmp0 <- FindNeighbors(tmp0, reduction = "harmony", dims = 1:dims)#harmony

res<-0.1 #0.1 resolution
tmp0 <- FindClusters(tmp0,resolution =res)
current.cluster.ids <- levels(Idents(tmp0))
new.cluster.ids <- as.numeric(current.cluster.ids)+1
Idents(tmp0) <- plyr::mapvalues(x = Idents(tmp0), from = current.cluster.ids, to = new.cluster.ids)
tmp0@meta.data$clusters<-Idents(tmp0)[rownames(tmp0@meta.data)]
table(tmp0@meta.data$clusters)

# umap
p1 <- DimPlot(tmp0, reduction = "umap", label = F,group.by = "Sample",label.size = 4,pt.size = 0)
p2 <- DimPlot(tmp0, reduction = "umap", label = T,group.by = "clusters",label.size = 4,pt.size = 0)
p<-plot_grid(p1,p2,ncol = 2)
pdf(paste0(prefix,'_sample_clusters_umap.pdf'),w=12,h=5)
print(p)
dev.off()


# rename marker gene 
genelist<-'PTPRC,CD3E,CD3D,CD3G,CD4,CD8A,CD8B,NCAM1,NKG7,GNLY,KLRF1,MKI67,CD79A,CD79B,CD19,MS4A1,IGHG1,MZB1,IGLC2,CSF1R,CD68,CD163,CD14,ITGAX,CD1C,S100A8,S100A9,FCGR3B,CSF3R,MS4A2,TPSAB1,TPSB2,CPA3,EPCAM,KRT19,PECAM1,VWF,FN1,COL1A1'
genes<-unique(c(unlist(strsplit(genelist,split = ",",fixed=T))))
names(genes)<-c(rep('T&NK',12),rep('B cell',7),rep('Mye',6),rep('Neu',4),rep('Mast',4),rep('Epi',2),rep('Endo',2),rep('Fib',2))
Idents(tmp0)<-'clusters'
tmp0 <- ScaleData(tmp0, features = genes, verbose = FALSE)
pdf('genes_clusters_dotplot.pdf',w=10,h=6,onefile=F)
DotPlot(tmp0,features=genes,cols = c("lightgrey", "red"))+theme(axis.text.x = element_text(size=10,angle =60, hjust = 1),legend.text = element_text(face="plain", color="black",size=8),legend.title =element_text(face="plain", color="black",size=10))+xlab('')+ylab('')
dev.off()



#3.find marker
tmp<-tmp0;prefix<-'BC'
Idents(tmp)<-'clusters'
cluster.averages <- AverageExpression(object = tmp, assays ='RNA',return.seurat = F)
write.csv(cluster.averages$RNA,file=paste(prefix,'clusters_averageExpression.csv',sep='_'),quote=F,row.names=T)
average.genes<-cluster.averages$RNA
all.markers<-FindAllMarkers(tmp, only.pos = T, min.pct =0.1, logfc.threshold = 0.25,test.use='wilcox')
markers<-subset(all.markers,p_val_adj<0.05)
da<-arrange(markers,cluster,desc(avg_log2FC))
write.csv(da,file=paste0(prefix,'_all.markers_FDR0.05.csv'),quote=F,row.names=T)
topn<-10
top_gene <- da %>% group_by(cluster) %>% top_n(topn,avg_log2FC)
filename<-paste0(prefix,'_top10_markers_genes_pheatmap.pdf')
genes<-unique(as.vector(top_gene$gene))
data<-average.genes
da0<-log2(data[genes,]+1)
n1<-length(unique(tmp$clusters))
pheatmap::pheatmap(da0,filename=filename,width=4+0.3*n1,height=3+0.8*n1,cluster_rows=F,cluster_cols=F,display_numbers = F,number_format = "%.0f",fontsize = 8,fontsize_col = 10,border_color=NA,angle_col = "45",scale='row')

#top5
da %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
tmp <- ScaleData(tmp, features = genes, verbose = FALSE)
pdf(paste0(prefix,'_top5_markers_genes_DotPlot.pdf'), w=10+n1,h=6+0.2*n1)
DotPlot(tmp,group.by='clusters', features=unique(top5$gene),cols = c("lightgrey", "red"))+theme(axis.text.x = element_text(size=10,angle =60, hjust = 1),legend.text = element_text(face="plain", color="black",size=8),legend.title =element_text(face="plain", color="black",size=10))+xlab('')+ylab('')
dev.off()


##---------------------------------------
#4. We subset TNK,B,Mye,Neu,Mast,Endo,Epi and Fib cells,and then recluster each celltype with harmony debatch effect method similared with major cluster in the front part except with parameters below...
prefix<-'TNK';nfeatures<-2000;npc<-30;dims<-20;res<-1.2
prefix<-'Mye';nfeatures<-2000;npc<-30;dims<-20;res<-0.5
prefix<-'Neu';nfeatures<-2000;npc<-30;dims<-20;res<-0.1
prefix<-'Mast';nfeatures<-1000;npc<-20;dims<-10;res<-0.1
prefix<-'Epi';nfeatures<-2000;npc<-30;dims<-20;res<-0.1
prefix<-'Endo';nfeatures<-1000;npc<-20;dims<-10;res<-0.1
prefix<-'Fib';nfeatures<-1000;npc<-20;dims<-10;res<-0.1

