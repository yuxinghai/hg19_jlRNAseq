#!/usr/bin/Rscript
setwd("/data3/zhoulab/yuxinghai/jlRNAseq/02_mapping")
directory <- '/data3/zhoulab/yuxinghai/jlRNAseq/02_mapping'
sampleFiles <- grep("_NEAT1_V2_featureCounts.txt",list.files(directory),value=TRUE)
for (i in sampleFiles) {
  file1 <-read.table(i)
  value <-file1[2,][2]-file1[1,][2]
  v1 <-c("NEAT1_v1",as.character(value))
  v2 <-c("NEAT1_V2",as.character(file1[1,][2]))
  NEAT1 <-rbind(v1,v2)
  write.table(NEAT1,paste(directory,paste("re",i,sep = "_"),sep = "/"),quote = F,col.names = F,row.names = F,sep = "\t")
}

