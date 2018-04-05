library(ggplot2)
library(RColorBrewer)
library(plyr)
library(Cairo)
library(ggthemes)
go_plot <-function(h4_cc_plot,name){
  h4_cc_plot <-h4_cc_plot[seq(1,15),]
  cc_name <-gsub("GO:\\d+~","",h4_cc_plot$Term,perl=T)
  cc_name <-gsub(",*","",cc_name,perl=T)
  # order 
  h4_cc_plot$cc_name <-factor(cc_name,levels = cc_name[order(-h4_cc_plot$PValue)])
  
  plot_ano <-theme(legend.position="right",
                   panel.grid=element_blank(),
                   legend.title = element_blank(),
                   legend.text= element_text(face="bold", color="black",family = "Times", size=8),
                   plot.title = element_text(hjust = 0.5),
                   axis.text.x = element_text(face="bold", color="black", size=12),
                   axis.text.y = element_text(face="bold",  color="black", size=12),
                   axis.title.x = element_text(face="bold", color="black", size=12),
                   axis.title.y = element_text(face="bold",color="black", size=12))
  
  p <-ggplot(h4_cc_plot)+
    theme_bw(base_size = 12, base_family = "Times") +
    geom_hline(yintercept=c(5,10), colour="grey", linetype="dotted")+
    geom_bar(aes(x=cc_name,y=-log10(PValue)),fill="#FF0000",width=0.3,position = "dodge", stat="identity")+
    geom_text(aes(label=Count,x=cc_name,y=-log10(PValue)),vjust=0.5,hjust=0,size = 5)+
    labs(title = "",y="-log10(PValue)",x="Term")+
    plot_ano+
    coord_flip()
  name=gsub("\\.tsv","",name)
  o_name <-file.path(paste(name,".png",sep=""))
  Cairo(file=o_name, 
        type="png",
        units="in",
        bg="white",
        width=8, 
        height=5, 
        pointsize=12, 
        dpi=300)
  print(p)
  dev.off()
  
}

