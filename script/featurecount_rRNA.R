#!/usr/bin/Rscript
library(Rsubread)
rRNA <- '/data3/zhoulab/yuxinghai/anotiation/grch38.rmsk.rna_diy.gtf'
name <-c("F1")
DESdir <-'/data3/zhoulab/yuxinghai/jlRNAseq/02_mapping'
for (i in name) {
  bsname <-paste(i,"tran_map.sam",sep = ".")
  filename <-paste(DESdir,i,bsname,sep = "/")
  fc_SE <- featureCounts(filename,annot.ext=rRNA,isGTFAnnotationFile=T,isPairedEnd=TRUE,strandSpecific=2,nthreads=6,countMultiMappingReads=T)
  write.table(fc_SE$counts,paste(DESdir,paste(i,"rRNA_featureCounts.txt",sep = "_"),sep = "/"),quote = F,col.names = F,sep = "\t")
  write.table(fc_SE$stat,paste(DESdir,paste(i,"rRNA_featureStat.log",sep = "_"),sep = "/"),quote = F,col.names = F,sep = "\t")
  
}

