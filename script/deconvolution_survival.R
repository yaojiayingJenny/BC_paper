#1.download TCGA-BLCA data set with GDCRNATools tools and creat tmp matrix
library(GDCRNATools)
library(DT)
gdcRNADownload(project.id = 'TCGA-BLCA',data.type  = 'RNAseq', directory='TCGA-BLCA/RNAseq')


#2.creat TCGA count matrix and BC_Minor_tumor.rds file
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-BLCA',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)
### filter sample
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
##head 
datatable(as.data.frame(metaMatrix.RNA[1:6,]), extensions = 'Scroller',
          options = list(scrollX = TRUE, deferRender = TRUE, scroller = TRUE))

write.table(metaMatrix.RNA,file='TCGA-BLCA.metaMatrix.RNA.csv',quote=F,sep='\t',row.names=F)


#dataset<-read.csv('TCGA-BLCA.metaMatrix.RNA.csv',header=T,sep='\t')
dataset<-metaMatrix.RNA
inpath<-'TCGA-BLCA/RNAseq/' #the path of downloaded RNAseq files
data<-subset(dataset,sample_type=='PrimaryTumor')
da<-read.table(filename<-paste0(inpath,data[1,2],'/',data[1,1]),sep='\t',skip=6)
tmp<-data.frame(gene_id=da$V1,gene_name=da$V2,gene_type=da$V3)
for (i in seq(1:nrow(data))){
filename<-paste0(inpath,data[i,2],'/',data[i,1])
print(i)
da<-read.table(filename,sep='\t',skip=6)
colnames(da)<-c('gene_id','gene_name','gene_type','unstranded','stranded_first','stranded_second','tpm_unstranded','fpkm_unstranded','fpkm_uq_unstranded')
#d<-as.vector(da[,'tpm_unstranded']) #TCGA-BLCA_406_log2TPM_genesymbol.csv
d<-as.vector(da[,'unstranded'])
tmp<-cbind(tmp,d)
}
colnames(tmp)<-c('gene_id','gene_name','gene_type',as.vector(data$sample))
norm<-tmp
norm1<-norm[!duplicated(norm$gene_name),]
rownames(norm1)<-norm1$gene_name
norm2<-norm1[,as.vector(data$sample)]
write.csv(norm2,file='TCGA-BLCA_406_gene_counts.csv',quote=F,row.names=T)

#BC_Minor_tumor.rds
tpm_norm<-norm2
bk.dat<-data.frame(t(tpm_norm))
saveRDS(bk.dat,file='bk.dat_TCGA-BLCA_406_log2TPM_genesymbol.rds')
rds0<-readRDS('Fig1/BC_celltype.rds') #major celltype rds
#subset with Tumor samples
rds<-subset(rds0,Group=='Tumor')
saveRDS(rds,file='Fig1/BC_Minor_tumor.rds')


#3.deconvolution with BayesPrism and HR survival analysis
#https://github.com/yaojiayingJenny/deconvolution_survival 
mkdir result
make -f mk_deconvolution_survival rds=BC_Minor_tumor.rds outdir=result tumor=TCGA-BLCA count=TCGA-BLCA_406_gene_counts.csv meta=TCGA-BLCA.metaMatrix.RNA.csv idents=Minor deconvolution_survival

