directory <-"/data2/zhoulab/yuxinghai/hg19_jlRNAseq/results/03_featurecount"
# directory <-"~/mount_dir/rna1_data2/hg19_jlRNAseq/results/03_featurecount"
setwd(directory)
#distance in heatmap(all samples)
library(DESeq2)
library(Cairo)
library(corrplot)
library("RColorBrewer")
#install.packages("corrplot")
sampleFiles <- grep("_featureCounts.txt",list.files(directory),value=TRUE)
#head(sampleFiles)
sampleCondition <-c("ipsf","ipsf","control","control","inono","inono","iNT","iNT","ipspc1"
                    ,"ipspc1","iV2","iV2")

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition) 
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('control','ipsf',"inono","iNT"
                                                                          ,"ipspc1","iV2"))
dds<-DESeq(ddsHTSeq)
colnames(dds) <-c("F1","F2","G1","G2","NO1","NO2","NT1","NT2","PS1","PS2","V2_1","V2_2")
res<-results(dds)
sampleCor <- cor(assay(dds))
#colnames(sampleDistMatrix) <- NULL   

Cairo(file="figure/distance_PNG_300_lengend_dpi.png", 
      type="png",
      units="in",
      bg="white",
      width=8, 
      height=8, 
      pointsize=12, 
      dpi=300)

p <-corrplot(sampleCor, method="shade",addshade= "positive" ,shade.col=NA, tl.col="black", tl.srt=45,
             addCoef.col="black", addcolorlabel="no", order="alphabet",type = "upper",cl.lim=(0:1))

print(p)
dev.off()
