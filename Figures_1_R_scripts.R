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

rm(list = ls())

############set working file folder############
setwd(dir = "E:/immunogenemoci landsapes in subsites of HNSCC/manuscript/1.26/Figures_from_R/")


##################Figure 1 HPV RNA and P16################

load(file = 'E:/immunogenemoci landsapes in subsites of HNSCC/immunity score and genes/clin.data.HPV.HNSC.CESC.RData')  ####load clinical data


table(clin$Smoking,clin$HPV_status_Subsite.2)

clin$Smoking=recode(clin$Smoking,'Current reformed smoker for < or = 15 years'='Current reformed smoker' 
                    
                    ,'Current reformed smoker for > 15 years'='Current reformed smoker'                        
                    
                    ,'Current Reformed Smoker, Duration Not Specified'='Current reformed smoker' )

HPV.load=readxl::read_xlsx(path = 'E:/immunogenemoci landsapes in subsites of HNSCC/immunity score and genes/raw data/cell report/table1 HPV status.xlsx',
                           sheet = 'Table S1P',col_names = T)%>%select(
                             ParticipantBarcode =1,
                             BC.RNA.HPV_count =5,
                             BC.RNA.HPV_rpm=6,
                             BC.exome.HPV_count=11,
                             BC.exome.HPV_rpm=12
                           )           ####load HPV RNA expression and exome 

HPV.load=HPV.load[-c(1:2),]

HPV.load[,2:5]=apply(HPV.load[,2:5],2,as.numeric)


###########################P16 information######

p16_hpv_ish_TCGA=read.csv('F:/TCGA/TCGA trainning/HNSCC/clinical_trimmed_data_527cancer.csv',header = T
)%>%select(
  'Patient.ID',
  'Hpv.status.ish',
  'Hpv.status.p16'
)

colnames(p16_hpv_ish_TCGA)=c('ParticipantBarcode','HPV_ISH','P16_IHC')

p16_hpv_ish_TCGA$P16_IHC[p16_hpv_ish_TCGA$P16_IHC==""]="No_test"

p16_hpv_ish_TCGA$HPV_ISH[p16_hpv_ish_TCGA$HPV_ISH==""]="No_test"

##################E6 and E7 expression from sulin paper#####
{
sulin_E1234567=readxl::read_xlsx(path = 'F:/TCGA/TCGA trainning/HNSCC/sulin_HPVE1234567 expression HNSC.xlsx',
                                 col_names = F)

sulin_E1234567=sulin_E1234567[-1,]
rn=sulin_E1234567$...1[2:nrow(sulin_E1234567)]
cn=sulin_E1234567[1,2:ncol(sulin_E1234567)]
cn=toupper(cn)
sulin_E1234567=sulin_E1234567[-1,-1]

colnames(sulin_E1234567)=cn

sulin_E1234567=apply(sulin_E1234567, 2, as.numeric)

rownames(sulin_E1234567)=rn

table(colSums(sulin_E1234567)>0)

E_df=sulin_E1234567[,colSums(sulin_E1234567)>0]%>%t()%>%as.data.frame()

#E_df$HPV_sums_E=rowSums(E_df)

E_df$ParticipantBarcode=rownames(E_df)

E_df1=E_df%>%melt(id.vars="ParticipantBarcode")

E_df1$variable=str_split(E_df1$variable,pattern = "_",simplify = T)[,2]

E_df1=E_df1%>%filter(value!=0)%>%dcast(ParticipantBarcode~variable, value.var = "value", fun.aggregate = sum)
}
########################

clin.hpv=merge(clin,HPV.load)

meta1=clin.hpv[,c(1,9,10,34:40)]%>%merge(p16_hpv_ish_TCGA,.)

meta1$P16_IHC=factor(meta1$P16_IHC,levels = c("Negative","Positive","No_test"))

meta2=E_df1%>%merge(.,clin.hpv[,c(1,9,10,34:40)])

meta1$Subsite.2=factor(meta1$Subsite.2,levels = c("OP","Non-OP"))

meta1$Subsite.1=factor(meta1$Subsite.1,levels = c("Oropharynx","Oral.cavity","Larynx","Hypopharynx"))

########################P16 is ok for OPSCC but not in Non_OPSCC###############################

tmp=meta1%>%select(ParticipantBarcode,Subsite.2,Subsite.1,P16_IHC,BC.RNA.HPV_rpm)


tmp2=clin.hpv%>%filter(HPV_status=='positive')


tt=tapply(tmp2$BC.RNA.HPV_rpm,tmp2$Subsite.1,summary)


tt1=bind_rows(tt)
rownames(tt1)=names(tt)


tt=tapply(tmp2$BC.RNA.HPV_rpm,tmp2$Subsite.2,summary)

tt1=bind_rows(tt)
rownames(tt1)=names(tt)



tt=tapply(clin.hpv$OS.time,clin.hpv$Subsite.2,summary)

tt1=bind_rows(tt)
rownames(tt1)=names(tt)



tt=tapply(clin.hpv$OS.time,clin.hpv$Site,summary)

tt1=bind_rows(tt)
rownames(tt1)=names(tt)






tmp$BC.RNA.HPV_rpm=log10(tmp$BC.RNA.HPV_rpm+1)

non_OP_HNSC=tmp%>%filter(Subsite.2=="Non-OP")

colnames(non_OP_HNSC)[2]='Subsite'

non_OP_HNSC=non_OP_HNSC[,-3]

HNSC_p16=tmp[,-2]

colnames(HNSC_p16)[2]='Subsite'

tmp1=rbind(HNSC_p16,non_OP_HNSC)

tmp1$Subsite=factor(tmp1$Subsite,levels = c("Oropharynx",
                                            "Non-OP",
                                            "Oral.cavity",
                                            "Larynx",
                                            "Hypopharynx"))

for_label=tmp1%>%filter(P16_IHC=="Positive"&Subsite=="Non-OP"&BC.RNA.HPV_rpm>1)


F1=tmp1%>%select(P16_IHC,Subsite,BC.RNA.HPV_rpm)%>%
  ggboxplot(., x = "P16_IHC", y ="BC.RNA.HPV_rpm" ,
            facet.by = 'Subsite',
            add = "jitter",
            add.params = list(size=0.5),
            scales='free_y',
            color = "P16_IHC",
            ylab = "HPV RNA load Log10(Reads per Million+1)",
            xlab = "HPV status determined by P16 IHC in HNSC",
            order = c("Negative","Positive","No_test"),
            ncol=4,
            lwd=0.2,
            width = 0.7,
            panel.labs = list(Subsite = c("A. OPSCC","B. non-OP HNSC","C. OCSCC","D. LSCC","E. HPSCC")), #####facet legend text
            panel.labs.font = list(face = 'bold', #color = 'white', 
                                   size = 10, angle = NULL,family='sans',hjust=0), ####facet legend text font
            palette = "aaas",#fill = 'HPV_type',
            legend="none"
  )+
  #scale_y_continuous(expand = expansion(mult = c(0.05, 1)))+
  #yscale("log10", .format = TRUE)+
  annotation_logticks(sides = "l")+
  rotate_x_text(angle = 45)+
  theme(strip.background =element_rect(fill = "transparent",  #####facet background
                                       colour = NA))+
  stat_compare_means(comparisons = list(c("Negative",'Positive')),label = "p.format",label.y.npc = 'top',size=3)+
  font("xy", size=10,family = 'sans',face = "bold")+
  font("axis.text", size=8,family = 'sans')+
  font("axis.title", size=10,family = 'sans')+
  font("legend.title", size=10,family = 'sans')+
  font("legend.text", size=10,family = 'sans')+
  scale_x_discrete(labels=c('P16 negative','P16 positive','Unkown'))+ 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))+  ####margin between top and bottom main and xy axis
  font("title", size=10,family = 'sans')+
 # geom_hline(yintercept = log10(1.2),linetype="dashed", color="red")+
  ggsave(filename = "./Figure 1/HPV_RNA_load_in_P16_more_subsites.png",
         width = 10,height = 5,units = "in")


#####################################
meta3=tmp1%>%left_join(.,E_df1)

var1=c("E1","E2","E4","E5","E6","E7","L1","L2")

meta3[,var1][is.na(meta3[,var1])]=0

meta3[,var1]=log10(meta3[,var1]+1)

meta3=meta3%>%filter(!is.na(Subsite))

F2=meta3%>%select(Subsite,P16_IHC,E6,E7)%>%filter(Subsite%in%c("Oropharynx","Non-OP"))%>%
  filter(P16_IHC!="No_test")%>%
  melt(id.vars=c("Subsite","P16_IHC"))%>%
  ggviolin(.,x="variable",y="value",
           color = "P16_IHC",facet.by = "Subsite"
           ,palette = "aaas"
           ,xlab = "HPV related biomarkers"
           ,ylab = "RNA expression log10(Read counts+1)"
           ,legend="none"
           ,panel.labs = list(Subsite = c("F. OPSCC","G. non-OP HNSC","H. OCSCC","I. LSCC","J. HPSCC")), #####facet legend text
           panel.labs.font = list(face = 'bold', #color = 'white', 
                                  size = 10, angle = NULL,family='sans',hjust=0) ####facet legend 
           ,ncol=1
           ,size = 0.5
           ,width = 1
           , add="jitter",
           add.params = list(size=0.3),
           
  )+stat_compare_means(aes(group = P16_IHC),
                       label = "p.format",label.y.npc = 'top',size=3)+
  theme(strip.background =element_rect(fill = "transparent",  #####facet background
                                       colour = NA))+
  font("xy", size=10,family = 'sans',face = "bold")+
  font("axis.text", size=8,family = 'sans')+
  font("axis.title", size=10,family = 'sans')+
  font("legend.title", size=10,family = 'sans')+
  font("legend.text", size=10,family = 'sans')+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))+  ####margin between top and bottom main and xy axis
  font("title", size=10,family = 'sans')+
  ggsave(filename = "./Figure 1/F1_2_Es_in_P16_IHC.png",width = 4,height = 3,units = "in")



F2_1=meta3%>%select(Subsite,P16_IHC,E6,E7)%>%filter(Subsite!="Hypopharynx")%>%
  filter(P16_IHC!="No_test")%>%
  melt(id.vars=c("Subsite","P16_IHC"))%>%
  ggviolin(.,x="variable",y="value",
           color = "P16_IHC",facet.by = "Subsite"
           ,palette = "aaas"
           ,xlab = "HPV related biomarkers"
           ,ylab = "RNA expression log10(Read counts+1)"
           ,legend="none"
           ,panel.labs = list(Subsite = c("F. OPSCC","G. non-OP HNSC","H. OCSCC","I. LSCC","J. HPSCC")), #####facet legend text
           panel.labs.font = list(face = 'bold', #color = 'white', 
                                  size = 10, angle = NULL,family='sans',hjust=0) ####facet legend 
           ,ncol=1
           ,size = 0.5
           ,width = 1
           , add="jitter",
           add.params = list(size=0.3),
           
  )+stat_compare_means(aes(group = P16_IHC),
                       label = "p.format",label.y.npc = 'top',size=3)+
  theme(strip.background =element_rect(fill = "transparent",  #####facet background
                                       colour = NA))+
  font("xy", size=10,family = 'sans',face = "bold")+
  font("axis.text", size=8,family = 'sans')+
  font("axis.title", size=10,family = 'sans')+
  font("legend.title", size=10,family = 'sans')+
  font("legend.text", size=10,family = 'sans')+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))+  ####margin between top and bottom main and xy axis
  font("title", size=10,family = 'sans')+
  ggsave(filename = "./Figure 1/F1_2_Es_in_P16_IHC_1.png",width = 4,height = 6,units = "in")



s_F1=meta3%>%select(Subsite,P16_IHC,var1)%>%
  filter(P16_IHC!="No_test")%>%
  melt(id.vars=c("Subsite","P16_IHC"))%>%
  ggviolin(.,x="variable",y="value",
           color = "P16_IHC",facet.by = "Subsite"
           ,palette = "aaas"
           ,xlab = "HPV related biomarkers"
           ,ylab = "RNA expression log10(Read counts+1)"
           ,legend="none"
           ,panel.labs = list(Subsite = c("A. OPSCC","B. non-OP HNSC","C. OCSCC","D. LSCC","E. HPSCC")), #####facet legend text
           panel.labs.font = list(face = 'bold', #color = 'white', 
                                  size = 10, angle = NULL,family='sans',hjust=0) ####facet legend 
           ,ncol=1
           ,size = 0.5
           ,width = 1
           , add="jitter",
           add.params = list(size=0.3),
           
  )+stat_compare_means(aes(group = P16_IHC),
                       label = "p.format",label.y.npc = 'top',size=3)+
  theme(strip.background =element_rect(fill = "transparent",  #####facet background
                                       colour = NA))+
  font("xy", size=10,family = 'sans',face = "bold")+
  font("axis.text", size=8,family = 'sans')+
  font("axis.title", size=10,family = 'sans')+
  font("legend.title", size=10,family = 'sans')+
  font("legend.text", size=10,family = 'sans')+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))+  ####margin between top and bottom main and xy axis
  font("title", size=10,family = 'sans')+
  ggsave(filename = "./Figure 1/Es_in_P16_IHC.pdf",width = 8,height = 8,units = "in")


ggarrange(ggarrange(F1,ncol = 1),ggarrange(F2,NULL,ncol = 2,widths = c(1.7,3)),
          nrow = 2,
          heights = c(1,1))+
  ggsave(filename ="./Figure 1/new_figure1.pdf",width = 10,height = 7,units = "in" )



############################################### ROC curves######################

library(pROC)

meta4=meta1%>%filter(P16_IHC%in%c("Negative","Positive"))

meta4$p16=ifelse(meta4$P16_IHC=="Negative",0,1)

meta4$HPV_status=factor(meta4$HPV_status,levels = c("negative","positive"))

meta4$HPV_status_1=ifelse(meta4$HPV_status=="negative",0,1)

meta5=split(meta4,meta4$Subsite.2)
roc_i<- roc(meta5$`Non-OP`,
            response = "HPV_status_1",
            predictor = 'p16',
            legacy.axes=TRUE,
            ci=TRUE, 
            boot.n=2000, 
            ci.alpha=0.9, 
            stratified=TRUE,
            plot=TRUE, 
            auc.polygon=TRUE, 
            max.auc.polygon=TRUE, 
            grid=TRUE,
            print.auc=TRUE, 
            print.thres='best', 
            print.thres.best.method='y',
            print.thres.adj=c(-0.05, 1.25),#,c(-0.05, 2.0),c(-0.05, 5.0)),
            print.thres.pattern="Cut-off: %.3f \n\nSp: %.3f \nSe: %.3f",
            print.thres.pattern.cex = 1.5,
            print.auc.cex = 1.5,
            print.auc.y=0.2,
            print.auc.x=0.8,
            cex.axis=1.5,
            cex.lab=1.5,
            print.thres.pch=16,
            print.thres.cex=2.0,
            cex.main=1.5)


roc_i2 <- roc(predictor = meta5$OP$p16,
              response = meta5$OP$HPV_status,
              legacy.axes=TRUE,
              ci=TRUE, 
              boot.n=2000, 
              ci.alpha=0.9, 
              stratified=TRUE,
              plot=TRUE, 
              auc.polygon=TRUE, 
              max.auc.polygon=TRUE, 
              grid=TRUE,
              print.auc=TRUE, 
              print.thres='best', 
              print.thres.best.method='y',
              print.thres.adj=c(-0.05, 1.25),#,c(-0.05, 2.0),c(-0.05, 5.0)),
              print.thres.pattern="Cut-off: %.3f \n\nSp: %.3f \nSe: %.3f",
              print.thres.pattern.cex = 1.5,
              print.auc.cex = 1.5,
              print.auc.y=0.2,
              print.auc.x=0.8,
              cex.axis=1.5,
              cex.lab=1.5,
              print.thres.pch=16,
              print.thres.cex=2.0,
              cex.main=1.5)

roclist <- list("P16 in non-OPSCC ROC" = roc_i,
                "P16 test in OPSCC ROC" = roc_i2
)


F3=ggroc(roclist, 
         aes = c("linetype", "color"),
         size=1.2,
         #alpha=0.5,
         legacy.axes = TRUE) +
  geom_abline(size=1,alpha=0.3,linetype="longdash") +
  theme_pubr() +
  ggtitle("E. ROC curves for P16 test predicting HPV status in OPSCC and non-OPSCC") +
  labs(x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)",
       linetype = "susbites"
  )+theme(
    legend.position = "none",
    legend.title = element_blank(),
    title = element_text(size=10,face = "bold",family = "sans"),
    legend.text = element_blank(),
    axis.title = element_text(size=8,face = "bold",family = "sans"),
    axis.text = element_text(size=8,family = "sans"),
    axis.line = element_line(size=1)
  )+scale_color_manual(values=c('#999999','#E69F00'))+
  scale_linetype_manual(values=c("solid", "solid"))+
  geom_label(
    aes(label = "OPSCC (AUC: 0.922, Specificity: 87.5%, Sensitivity: 96.9%)",
        x = 0.105,y=0.72),
    color="white",
    fill="#E69F00",
    nudge_y = 0.3,
    #nudge_x = -0.4,
    family = "sans",
    vjust=0,
    hjust=0.1,
    size=2.5,
    fontface = "bold"
    
  )+
  geom_label(
    aes(label = "non-OPSCC (AUC: 0.629, Specificity: 92.5%, Sensitivity: 33.3%)",
        x = 0.105,y=0.9),
    color="white",
    fill="#999999",
    nudge_y = 0.3,
    #nudge_x = -0.4,
    vjust=0,
    hjust=0.1,
    size=2.5,
    fontface = "bold"
  )+geom_text(aes(x=.37, y=1.129, label="VS.",
                  size=0.5,fontface=3),
              color="red")+
  ggsave(filename = "./Figure 1/ROC1.png")


