base <-"~/mount_dir/rna1_data2/hg19_jlRNAseq/results/03_featurecount/david"
setwd(base)
sampleFiles <- grep("tsv",list.files(base),value=TRUE)
# Import go_plot function
# df,name,out_dir
source("~/mount_dir/rna1_data2/hg19_jlRNAseq/script/david_plot.R")
for (samps in sampleFiles) {
  psf_bp <-read.delim(paste(base,samps,sep = "/"),header = T,sep = "\t")
  go_plot(psf_bp,samps)
}


