library(Seurat)
library(sctransform)
library(dplyr)
library(monocle3)
#library(monocle)
library(ggplot2)
library(cowplot)
library(configr)


cds2<-readRDS('Mono_Mac_tumor_immune_combined_monocle3.rds')
Mye_colors<-c('#D51F26','#208046','#2F2D66','#6E4B9E','#FBCB0A','#D24B27','#C06CAB','#D8A767','#89C75F','#90D5E4','#F37B7D','#9983BD','#8A9FD1')
# c('Mono1_CD14','Mono2_CD14','Mono_CD16','Mac_TREM2','Mac_CCL3','Mac_IL1B','Mac_CD163')
# c('cDC_CD1C','cDC_CD207','cDC_CLEC9A','DC_LAMP3')
prefix<-'Mono_Mac_tumor'
#plot_cells(lung, markers="GDF15")
monocle.object<-cds2
pdf(paste(prefix,'monocle3_celltype.pdf',sep='_'),w=5.5,h=4)
plot_cells(monocle.object,
           color_cells_by = "Minor",
           label_groups_by_cluster=FALSE,
           label_leaves=F,
           trajectory_graph_color = "black",
           label_roots = F,         
           cell_size = 0.5,rasterize =F,cell_stroke =NA,
           trajectory_graph_segment_size = 0.6,group_label_size  = 4,show_trajectory_graph = F,label_branch_points = F)+scale_colour_manual(values=Mye_colors,name = 'Mono_Mac')+theme(legend.position='right')
dev.off()

pdf(paste(prefix,'order_cells_pseudotime.pdf',sep='_'),w=5,h=4)
plot_cells(monocle.object, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,label_branch_points = FALSE,show_trajectory_graph = T,label_roots = F)
dev.off()
#+ scale_colour_continuous(low='#1D068D',high ='#F4EA26', name = 'Pseudotime',guide=guide_colorbar(reverse=F))


dataset<-read.table('Mono_Mac_tumor_Trajectory_DEG_genes.xls',header=T,sep='\t')
data<-arrange(dataset,desc(morans_test_statistic))

dim(pData(cds2))
genes<-data[1:10,]$gene_short_name
#Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
pdf(paste0(prefix,'_top10_Genes_Jitterplot.pdf'),width = 8, height = 6)
p <- plot_genes_in_pseudotime(monocle.object[genes,], color_cells_by='Minor', 
                             min_expr=0.5, ncol = 2,label_by_short_name = FALSE)+scale_colour_manual(values=Mye_colors,name = 'Mono_Mac')
print(p)
dev.off()
write.table(data,file='Mono_Mac_tumor_Trajectory_DEG_genes_sort.xls',quote=F,row.names=F,sep='\t')





##########DC-------------------------------------------------
cds2<-readRDS('DC3_tumor_immune_combined_monocle3.rds')
Mye_colors<-c('#D51F26','#208046','#2F2D66','#6E4B9E','#FBCB0A','#D24B27','#C06CAB','#D8A767','#89C75F','#90D5E4','#F37B7D','#9983BD','#8A9FD1')[8:11]
Mye_colors<-c('#D8A767','#89C75F','#F37B7D')
prefix<-'DC3_tumor'
#plot_cells(lung, markers="GDF15")
monocle.object<-cds2
pdf(paste(prefix,'monocle3_celltype.pdf',sep='_'),w=5.5,h=4)
plot_cells(monocle.object,
           color_cells_by = "Minor",
           label_groups_by_cluster=FALSE,
           label_leaves=F,
           trajectory_graph_color = "black",
           label_roots = F,         
           cell_size = 1,rasterize =F,cell_stroke =NA,
           trajectory_graph_segment_size = 0.6,group_label_size  = 4,show_trajectory_graph = F,label_branch_points = F)+scale_colour_manual(values=Mye_colors,name = 'DC')+theme(legend.position='right')
dev.off()

pdf(paste(prefix,'order_cells_pseudotime.pdf',sep='_'),w=5,h=4)
plot_cells(monocle.object, color_cells_by = "pseudotime",trajectory_graph_segment_size = 0.3, label_cell_groups = FALSE, label_leaves = FALSE,label_branch_points = FALSE,show_trajectory_graph = T,label_roots = F)
dev.off()
#+ scale_colour_continuous(low='#1D068D',high ='#F4EA26', name = 'Pseudotime',guide=guide_colorbar(reverse=F))


dataset<-read.table('DC3_tumor_Trajectory_DEG_genes.xls',header=T,sep='\t')
data<-arrange(dataset,desc(morans_test_statistic))

dim(pData(cds2))
genes<-data[1:10,]$gene_short_name
#Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
##差异基因表达趋势图
pdf(paste0(prefix,'_top10_Genes_Jitterplot.pdf'),width = 8, height = 6)
p <- plot_genes_in_pseudotime(monocle.object[genes,], color_cells_by='Minor', 
                             min_expr=0.5, ncol = 2,label_by_short_name = FALSE)+scale_colour_manual(values=Mye_colors,name = 'Mono_Mac')
print(p)
dev.off()

write.table(data,file='DC3_tumor_Trajectory_DEG_genes_sort.xls',quote=F,row.names=F,sep='\t')