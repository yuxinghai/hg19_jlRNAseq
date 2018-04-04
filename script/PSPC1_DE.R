base <-'/data2/zhoulab/yuxinghai/hg19_jlRNAseq/'
#base <-'~/mount_dir/rna1_data2/hg19_jlRNAseq/'

dir<-paste(base,'results/03_featurecount',sep="")
setwd(dir)

period <-"PSPC1"
genename_dir <- paste("gene_name",period,sep="/")
if (!file.exists(genename_dir)){
  dir.create(genename_dir)
}

figure_dir <- paste("figure",period,sep="/")
if (!file.exists(figure_dir)){
  dir.create(figure_dir)
}

DEgene_dir <- paste("DEgene",period,sep="/")
if (!file.exists(DEgene_dir)){
  dir.create(DEgene_dir)
}


# first 3 sample is control,later 3 is damage
library(DESeq2)
sampleFiles <- grep("_featureCounts.txt",list.files(dir),value=TRUE)
indx= c(3,4,seq(9,10))
h_12 <- sampleFiles[indx]
#head(sampleFiles)
sampleCondition <-c("control","control","KD","KD")

sampleTable <- data.frame(sampleName = h_12,
                          fileName = h_12,
                          condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       design= ~ condition,
                                       directory = dir)
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('control','KD'))
dds<-DESeq(ddsHTSeq)
colnames(dds) <-gsub("_featureCounts.txt","",h_12)
# no filter the gene count number 
res<-results(dds)
#write all expression gene name to a file
all_exp <-res[res$baseMean !=0,]
res<-res[order(res$padj),]
mcols(res, use.names=TRUE)
#head(res)

#add gene symbol and filter diff expressed gene
#up:log2FC>0.58, down:log2fc<-0.58
res <-res[complete.cases(res),] 
head(res)
res <-res[ res$padj <0.05,]
res$name <- gsub("\\.[0-9]+","",rownames(res))

new <-as.data.frame(res)
down_gene <- new[new$log2FoldChange<(-0.58),]
#sort by -log2fc
down_gene <-down_gene[with(down_gene,order(padj,-abs(log2FoldChange))),]

up_gene <- new[new$log2FoldChange>log2(1.5),]
up_gene <-up_gene[with(up_gene,order(padj,-abs(log2FoldChange))),]
up_and_down<-new[(new$log2FoldChange<(-0.58)) | (new$log2FoldChange>log2(1.5)), ]
#genename for go term

write.table(up_and_down$name,paste(genename_dir,"sig",sep="/"),sep = "\n",quote = F,col.names = F,row.names = F)
write.table(down_gene$name,paste(genename_dir,"down",sep="/"),sep = "\n",quote = F,col.names = F,row.names = F)
write.table(up_gene$name,paste(genename_dir,"up",sep="/"),sep = "\n",quote = F,col.names = F,row.names = F)


#volcan plot
source("../../script/volcan_plot.R")
cairo_plot(dds,paste(figure_dir,period,sep = "/"))

#DISTRIBUTION level of DE 

notsig <-dim(all_exp)[1]-dim(up_and_down)[1]
gene_num <-c(dim(up_gene)[1],notsig,dim(down_gene)[1])
expre_level <-c("Up","Not_change","Down")


plot <-data.frame(gene_num=gene_num,expre_level=expre_level)
library(ggplot2)
ggplot(plot)+geom_bar(aes(x=expre_level,y=gene_num),stat = "identity",
                      position = "dodge",width = 0.5)+labs(y="Gene number")+
  geom_text(aes(x=expre_level,y=gene_num,label=gene_num),vjust = -0.1 )
ggsave(paste(figure_dir,"DE_number.pdf",sep="/"))

#######################ano siggene####################

ano <- read.csv2("../../anno/hg19_david_anno.tsv",
                 sep="\t",header =T ,comment.char = "#")

down_gene <-merge(down_gene,ano,by.x="name",by.y="ID",all.x=T)
#head(down_gene)
up_gene <-merge(up_gene,ano,by.x="name",by.y="ID",all.x=T)
write.table(down_gene[,c(1,8,2,3,6,7,10,11,12,13,14,15)],paste(DEgene_dir,"ano_downgene.tsv",sep="/"),sep="\t",quote = F,row.names = F)
write.table(up_gene[,c(1,8,2,3,6,7,10,11,12,13,14,15)],paste(DEgene_dir,"ano_upgene.tsv",sep="/"),sep="\t",quote = F,row.names = F)

