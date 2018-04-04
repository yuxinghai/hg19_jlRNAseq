#!/usr/bin/Rscript
library(Rsubread)
GTF <- '/data2/zhoulab/yuxinghai/hg19_jlRNAseq/anno/gencode.v19.annotation.gtf'
name <-c("F1","F2","G1","G2","N01","N02","NT1","NT2","PS1","PS2","V1","V2")
DESdir <-'/data2/zhoulab/yuxinghai/hg19_jlRNAseq/results/02_mapping/sorted'
out_dir <-'/data2/zhoulab/yuxinghai/hg19_jlRNAseq/results/03_featurecount'
for (i in name) {
  bsname <-paste(i,"uniq_sort.bam",sep = "_")
  filename <-paste(DESdir,i,bsname,sep = "/")
  fc_SE <- featureCounts(filename,annot.ext=GTF,isGTFAnnotationFile=T,isPairedEnd=TRUE,strandSpecific
                         =2,nthreads=6)

  write.table(fc_SE$counts,paste(out_dir,paste(i,"featureCounts.txt",sep = "_"),sep = "/"),quote = F,col.names = F,sep = "\t")
  write.table(fc_SE$stat,paste(out_dir,paste(i,"featureStat.log",sep = "_"),sep = "/"),quote = F,col.names = F,sep = "\t")
  
}

