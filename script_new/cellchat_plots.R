library(CellChat)
library(ggplot2)
library(patchwork)
options(stringsAsFactors = FALSE)
library(cowplot) 
library(dplyr)
library(Seurat) #CombinePlots


###--------------------------------------------------------------------
#FigS6c
cellchat<-readRDS('Tumor_vs_Normal_cellchat_merge.rds')
pathways.show.all <- unique(c(cellchat@netP[[1]]$pathways,cellchat@netP[[2]]$pathways))
use1<-'Neutrophil';use2<-'Neutrophil';prefix<-'Tumor_vs_Normal'
gg1 <- rankNet(cellchat, mode = "comparison",sources.use = use1, targets.use = NULL,stacked = T, do.stat = TRUE,color.use=c('blue','red'))+ggtitle(paste0('sources=',use1)) ##group1,group2
gg2 <- rankNet(cellchat, mode = "comparison", sources.use = use1, targets.use = NULL,stacked = F, do.stat = TRUE,color.use=c('blue','red'))+ggtitle(paste0('sources=',use1))
gg3 <- rankNet(cellchat, mode = "comparison",sources.use = NULL, targets.use = use2,stacked = T, do.stat = TRUE,color.use=c('blue','red'))+ggtitle(paste0('targets=',use2)) ##group1,group2
gg4 <- rankNet(cellchat, mode = "comparison", sources.use = NULL, targets.use = use2,stacked = F, do.stat = TRUE,color.use=c('blue','red'))+ggtitle(paste0('targets=',use2))
gg<-gg1 + gg2+gg3 + gg4
pdf(paste0(prefix,'_source-',use1,'_target-',use2,'_comparison_pathway_rankNet.pdf'),w=10,h=0.1*length(pathways.show.all)+3)
print(gg)
dev.off()



pathways.show.all <- unique(c(cellchat@netP[[1]]$pathways,cellchat@netP[[2]]$pathways))
use1<-'Mast';use2<-'Mast';prefix<-'Tumor_vs_Normal'
gg1 <- rankNet(cellchat, mode = "comparison",sources.use = use1, targets.use = NULL,stacked = T, do.stat = TRUE,color.use=c('blue','red'))+ggtitle(paste0('sources=',use1)) ##group1,group2
gg2 <- rankNet(cellchat, mode = "comparison", sources.use = use1, targets.use = NULL,stacked = F, do.stat = TRUE,color.use=c('blue','red'))+ggtitle(paste0('sources=',use1))
gg3 <- rankNet(cellchat, mode = "comparison",sources.use = NULL, targets.use = use2,stacked = T, do.stat = TRUE,color.use=c('blue','red'))+ggtitle(paste0('targets=',use2)) ##group1,group2
gg4 <- rankNet(cellchat, mode = "comparison", sources.use = NULL, targets.use = use2,stacked = F, do.stat = TRUE,color.use=c('blue','red'))+ggtitle(paste0('targets=',use2))
gg<-gg1 + gg2+gg3 + gg4
pdf(paste0(prefix,'_source-',use1,'_target-',use2,'_comparison_pathway_rankNet.pdf'),w=10,h=0.1*length(pathways.show.all)+3)
print(gg)
dev.off()


###--------------------------------------------------------------------
#Fig6f
pathways.show <- c("ANNEXIN")
pathname<-'ANNEXIN'
colorlist<-c('#430D54','#FDE728')
pdf(paste0(pathname,'_netVisual_bubble_group_part_new.pdf'),w=15,h=4,onefile=T)
netVisual_bubble(cellchat, signaling = pathways.show,sources.use = c(1:9), targets.use = c(4,5),  comparison = c(1, 2), angle.x = 45,color.heatmap='viridis',n.colors=2,color.text=c('blue','red'),sort.by.target =  TRUE,thresh =1)
dev.off()