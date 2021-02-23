rm(list = ls())

############set working file folder############
setwd(dir = "E:/immunogenemoci landsapes in subsites of HNSCC/manuscript/1.26/Figures_from_R/Figure 3/")
library(data.table)
library(dplyr)
library(table1)
library(ggbeeswarm)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggstatsplot)
library(readxl)
library(survival)
library(survminer)
library(cowplot)
library(foreach)
library(showtext)
library(ComplexHeatmap)
library(gtsummary)

options(stringsAsFactors = F)


########################Fig3 ############################################

HPV.load=readxl::read_xlsx(path = 'E:/immunogenemoci landsapes in subsites of HNSCC/immunity score and genes/raw data/cell report/table1 HPV status.xlsx',
                           sheet = 'Table S1P',col_names = T)%>%select(
                             ParticipantBarcode =1,
                             BC.RNA.HPV_count =5,
                             BC.RNA.HPV_rpm=6,
                             BC.exome.HPV_count=11,
                             BC.exome.HPV_rpm=12
                           )
HPV.load=HPV.load[-c(1:2),]

HPV.load[,2:5]=apply(HPV.load[,2:5],2,as.numeric)

HPV=HPV.load%>%select(ParticipantBarcode,BC.RNA.HPV_rpm)

clin=merge(clin,HPV)


##########################immune transcript###################################################

imRNA.xy<-read.table(file = "E:/immunogenemoci landsapes in subsites of HNSCC/immunity score and genes/raw data/immune gene expression_HNSC_CESC.txt")

colnames(imRNA.xy)<-c('ParticipantBarcode','x','y')

imRNA.xy<-merge(imRNA.xy,clin)


#######
p<-
  ggplot(imRNA.xy,aes(x=y,y=x, colour=HPV_status#,shape=Subsite.2
  ))+
  geom_point(size=0.5,alpha=0.8)+
  scale_color_manual(values=c('blue','red'))+
  theme(title = element_text(size = 6),
        legend.position="top",
        legend.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=6,hjust = 0.5),
        panel.grid =element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        strip.background =element_rect(fill = "transparent",colour = NA),
        strip.text = element_text(size = 6),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )+ggsave(filename = "./immune transcript-tumormap.pdf",
           width = 40,height = 50,units = 'mm')


p

#########

p2<-p+facet_grid(~Subsite.2)

p2+ggsave(filename = "./facet-subsite-immune transcript-tumormap.pdf",
          width = 120,height = 50,units = 'mm')


p3=imRNA.xy%>%filter(Subsite.2=='Non-OP')%>%
  ggplot(aes(x=y,y=x, colour=HPV_status#,shape=Subsite.2
  ))+
  geom_point(size=0.5,alpha=0.8)+
  scale_color_manual(values=c('blue',"red"))+
  theme(title = element_text(size = 6),
        legend.position="top",
        legend.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=6,hjust = 0.5),
        panel.grid =element_blank(),
        axis.text = element_blank(),
        strip.text = element_text(size = 6),
        axis.ticks = element_blank(),
        strip.background =element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )+facet_grid(~Subsite.1)+
  ggsave(filename = "./facet-Non-OP-immune transcript-tumormap.pdf",
         width = 120,height = 50,units = 'mm')


library(patchwork)

(p2 / p3) +
  ggsave(filename = "./fig3-immune transcript-tumormap.1.pdf",
         width = 120,height = 115,units = 'mm')

####################################immune transcript supplemental figure 2###################################################
library(patchwork)
#########HPV subtype didn't matter########

imRNA.xy$HPV_count=(log(imRNA.xy$BC.RNA.HPV_rpm+1,base = 2))
imRNA.xy$HPV_subtype=imRNA.xy$`Major HPV type`
imRNA.xy$HPV_subtype[is.na(imRNA.xy$HPV_subtype)]='Negative'
imRNA.xy$HPV_subtype[!(imRNA.xy$HPV_subtype%in%c('HPV_16','HPV_18','Negative'))]='HPV_other'  

imRNA.xy$HPV_subtype=factor(imRNA.xy$HPV_subtype,
                            levels = c('Negative','HPV_16','HPV_other','HPV_18'))  
{    s2.a=ggplot(imRNA.xy,aes(x=y,y=x, colour=HPV_status#,shape=HPV_subtype
))+
    geom_point(size=0.5,alpha=0.7,aes(#color=HPV_count,
      shape=HPV_subtype))+
    #scale_color_gradient(low = 'blue',high = 'red')+
    scale_color_manual(values=c('blue','red'))+
    theme(title = element_text(size = 6),
          legend.position="top",
          legend.text = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size=6,hjust = 0.5),
          panel.grid =element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          strip.background =element_rect(fill = "transparent",colour = NA),
          strip.text = element_text(size = 6),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )+facet_grid(~Subsite.2)+ggsave(filename = "./Figure.S2.HPVsubtype1.pdf",
                                    width = 120,height = 50,units = 'mm')
  
  s2.b=imRNA.xy%>%filter(Subsite.2=='Non-OP')%>%
    ggplot(aes(x=y,y=x, colour=HPV_status#,shape=Subsite.2
    ))+
    geom_point(size=0.5,alpha=0.7,aes(#color=HPV_count,
      shape=HPV_subtype))+
    #scale_color_gradient(low = 'blue',high = 'red')+
    scale_color_manual(values=c('blue','red'))+
    theme(title = element_text(size = 6),
          legend.position="top",
          legend.text = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size=6,hjust = 0.5),
          strip.text = element_text(size = 6),
          panel.grid =element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.background =element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )+facet_grid(~Subsite.1)+
    ggsave(filename = "./Figure.S2.HPVsubtype2.pdf",
           width = 120,height = 50,units = 'mm')
  (s2.a/ s2.b) + #plot_layout(guides = 'collect')+
    ggsave(filename = "./Fig S2.HPV.subtype.pdf",
           width = 120,height = 115,units = 'mm')
}

##### Gender##
{    s2.a=ggplot(imRNA.xy,aes(x=y,y=x, colour=HPV_status#,shape=gender
))+
    geom_point(size=0.5,alpha=0.7,aes(#color=HPV_count,
      shape=gender))+
    #scale_color_gradient(low = 'blue',high = 'red')+
    scale_color_manual(values=c('blue','red'))+
    theme(title = element_text(size = 6),
          legend.position="top",
          legend.text = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size=6,hjust = 0.5),
          panel.grid =element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          strip.background =element_rect(fill = "transparent",colour = NA),
          strip.text = element_text(size = 6),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )+facet_grid(~Subsite.2)+ggsave(filename = "./Figure.S2.gender1.pdf",
                                    width = 120,height = 50,units = 'mm')
  
  s2.b=imRNA.xy%>%filter(Subsite.2=='Non-OP')%>%
    ggplot(aes(x=y,y=x, colour=HPV_status#,shape=Subsite.2
    ))+
    geom_point(size=0.5,alpha=0.7,aes(#color=HPV_count,
      shape=gender))+
    #scale_color_gradient(low = 'blue',high = 'red')+
    scale_color_manual(values=c('blue','red'))+
    theme(title = element_text(size = 6),
          legend.position="top",
          legend.text = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size=6,hjust = 0.5),
          strip.text = element_text(size = 6),
          panel.grid =element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.background =element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )+facet_grid(~Subsite.1)+
    ggsave(filename = "S2.gender2.pdf",
           width = 120,height = 50,units = 'mm')
  (s2.a/ s2.b) + #plot_layout(guides = 'collect')+
    ggsave(filename = "./Fig S2.gender.pdf",
           width = 120,height = 115,units = 'mm')
}
#######smoking####
##### Gender##
{    
  s2.a=ggplot(imRNA.xy,aes(x=y,y=x, colour=HPV_status#,shape=Smoking
  ))+
    geom_point(size=0.5,alpha=0.7,aes(#color=HPV_count,
      shape=Smoking))+
    #scale_color_gradient(low = 'blue',high = 'red')+
    scale_color_manual(values=c('blue','red'))+
    theme(title = element_text(size = 6),
          legend.position="top",
          legend.text = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size=6,hjust = 0.5),
          panel.grid =element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          strip.background =element_rect(fill = "transparent",colour = NA),
          strip.text = element_text(size = 6),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )+facet_grid(~Subsite.2)+ggsave(filename = "./Figure.S2.Smoking1.pdf",
                                    width = 200,height = 50,units = 'mm')
  
  s2.b=imRNA.xy%>%filter(Subsite.2=='Non-OP')%>%
    ggplot(aes(x=y,y=x, colour=HPV_status#,shape=Subsite.2
    ))+
    geom_point(size=0.5,alpha=0.7,aes(#color=HPV_count,
      shape=Smoking))+
    #scale_color_gradient(low = 'blue',high = 'red')+
    scale_color_manual(values=c('blue','red'))+
    theme(title = element_text(size = 6),
          legend.position="top",
          legend.text = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size=6,hjust = 0.5),
          strip.text = element_text(size = 6),
          panel.grid =element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.background =element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )+facet_grid(~Subsite.1)+
    ggsave(filename = "S2.Smoking2.pdf",
           width = 120,height = 50,units = 'mm')
  (s2.a/ s2.b) + #plot_layout(guides = 'collect')+
    ggsave(filename = "./Fig S2.Smoking.pdf",
           width = 120,height = 115,units = 'mm')
}


  
  

