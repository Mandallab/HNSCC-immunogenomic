rm(list = ls())

############set working file folder############
setwd(dir = "E:/immunogenemoci landsapes in subsites of HNSCC/manuscript/1.26/Figures_from_R/Figure 4/")
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

#########Fig 4########

load(file = 'F:/TCGA/pan-cancer/Pan-cancer/primary.solid.tumor.mRNA.expr.RData') ##mRNA expression data ps.data

ps.data$ParticipantBarcode<-row.names(ps.data)

mRNA<-ps.data%>%filter(ParticipantBarcode%in%clin$ParticipantBarcode) 

#rm(ps.data)

dat<-mRNA%>%select(-'ParticipantBarcode')%>%t()%>%apply(., 1, as.numeric)

rownames(dat)<-mRNA$ParticipantBarcode
rownames(mRNA)=mRNA$ParticipantBarcode
dat<-as.data.frame(t(dat))

res =immunedeconv::deconvolute(dat, "mcp_counter")

rownames(res)<-res$cell_type

res=res%>%as.data.frame()%>%select(-1)%>%t()%>%as.data.frame()

res$ParticipantBarcode=rownames(res)

save(res,file = 'MCP-Count.res.RData')

#############
immune.landscape<-read.csv(file = 'F:/TCGA/pan-cancer/Immunity score of TCGA cancer-Cell paper.csv')%>%select(-c(2,3,4))

colnames(immune.landscape)[1]<-"ParticipantBarcode"
immune.landscape=immune.landscape%>%mutate(TMB=Nonsilent.Mutation.Rate+Silent.Mutation.Rate,
                                           CNV.burden=Fraction.Altered+Number.of.Segments,
                                           Neoantigens=Indel.Neoantigens+SNV.Neoantigens)
colnames(res)<-c('T.cell',
                 'CD8.T.cell',
                 'Cytotoxicity.score',
                 'NK.cell',
                 'B.cell',
                 'Monocyte',
                 'Macrophage.Monocyte',
                 'Myeloid.dendritic.cell',
                 'Neutrophil',
                 'Endothelial.cell',
                 'Cancer.associated.fibroblast',
                 'ParticipantBarcode')


Immunosuppression.cytokine<-c('CXCL12','TGFB1','TGFB3','LGALS1')
T.cell.activation<-c('CXCL9','CXCL10','CXCL16','IFNG','IL15')
T.cell.survival<-c('CD70','CD27')
Regulatory.T.cell<-c('FOXP3','TNFRSF18')
MHC1<-c('HLA-A','HLA-B','HLA-C','HLA-C','HLA-E','HLA-F','HLA-G','B2M')
myeloid.chemotaxis<-c('CCL2')
TLS<-c('CCL19', 'CCL21','CXCL13', 'CCR7', 'CXCR5', 'SELL' , 'LAMP3',
       'CD79B','CD1D','CCR6','LAT','SKAP1','CETP','EIF1AY','RBP5','PTGDS')
Cytolytic.score=c("PRF1","GZMB")


gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

{res$Immunosuppression.cytokine<-apply(mRNA[,Immunosuppression.cytokine], 1, gm_mean)
  
  res$T.cell.activation<-apply(mRNA[,T.cell.activation], 1, gm_mean)
  
  res$T.cell.survival<-apply(mRNA[,T.cell.survival], 1, gm_mean)
  res$Regulatory.T.cell<-apply(mRNA[,Regulatory.T.cell], 1, gm_mean)
  res$MHC1<-apply(mRNA[,MHC1], 1, gm_mean)
  res$myeloid.chemotaxis<-mRNA[,myeloid.chemotaxis]
#  res$TLS<-apply(mRNA[,TLS.n1], 1, gm_mean)
  res$Cytolytic=apply(mRNA[,Cytolytic.score],1,gm_mean)
}


######################

meta1<-clin%>%
  select(ParticipantBarcode,Subsite.1,Subsite.2,contains('HPV_status'))%>%
  merge(.,res)%>%merge(.,immune.landscape)


######
var=c('B Cell Score','Naïve B Cell Proportion','BCR Diversity',
      'CD8+ T Cell', 
      'TCR Diversity',
      'Cytolytic activity score',
      'Cancer Associated Fibroblast Score',
      'Tumor Mutation Burden (per Mb)',
      'Copy Number Variation Burden')

df=meta1%>%select('Subsite.1',Subsite.2,'HPV_status',
                  'B Cell Score'='B.cell',
                  'Naive B Cell'='B.Cells.Naive',
                  'BCR diversity'='BCR.Shannon',
                  'CD8+ T Cell'='CD8.T.cell', 
                  'TCR diversity'='TCR.Shannon',
                  'Cytolytic activity score'='Cytolytic',
                  'Cancer Associated Fibroblast Score'='Cancer.associated.fibroblast',
                  'Tumor Mutation Burden (per Mb)'='TMB',
                  'Copy Number Variation Burden'='CNV.burden')

df$`B Cell Score`=log10(df$`B Cell Score`)
df$`CD8+ T Cell`=log10(df$`CD8+ T Cell`)
df$`Cytolytic activity score`=log10(df$`Cytolytic activity score`)  
df$`Cancer Associated Fibroblast Score`=log10(df$`Cancer Associated Fibroblast Score`)
df$`Copy Number Variation Burden`=log10(df$`Copy Number Variation Burden`)
df$`Tumor Mutation Burden (per Mb)`=log10(df$`Tumor Mutation Burden (per Mb)`)


df1=df%>%filter(Subsite.2=="Non-OP")

colnames(df1)[2]='Subsite'

df1=df1[,-1]

colnames(df)[1]='Subsite'

df=df[,-2]

df=rbind(df,df1)


reshape2::melt(df,id.vars =c('Subsite',
                             'HPV_status'),
               value.name = 'Immune.scores')%>%
  ggboxplot(., x = "Subsite", y ="Immune.scores" ,
            facet.by = 'variable',  
            panel.labs = list(variable = paste0(LETTERS[1:length(var)],'.  ',var)), #####facet legend text
            panel.labs.font = list(face = 'bold', #color = 'white', 
                                   size = 6, angle = NULL,family='sans',hjust=0), ####facet legend text font
            #panel.labs.background = list(color = "#cc80ff", fill = "#cc80ff", size = 0.5), #####facet legend background
            ncol=3, ###number of columns
            width = 0.7, ###box width
            #dodge=0.01,  ###distance between boxes
            scales='free_y', 
            color = 'HPV_status',
            shape="HPV_status",fill = '#F5F5F5', ###box background color
            outlier.shape = NA,     ####hide the outlier dots
            add.params=list(size=0,alpha=0.95), ###size and transparance 
            palette =c(
              "#63438c","#eb2ea5"), #fill = 'HPV_type', ###color for the group
            add = "jitter", ####the distribution of data 
            lwd=0.2, ##the box line width 
            # linetype = 'solidat', ###the box line type
            notch = F,
            order = c('Oropharynx',"Non-OP",'Oral.cavity','Larynx','Hypopharynx','Cervical'),
            font.label = list(size=6,color="black")
  )+rotate_x_text(angle = 45)+
  rotate_y_text(angle = 45)+
  scale_x_discrete(labels=c('OPSCC','non-OP HNSCC','OCSCC','LSCC','HPSCC','CESC'))+ 
  stat_compare_means(aes(group = HPV_status),label = "p.format",label.y.npc = 'top',size=1.5)+
  font("xy", size = 6,family = 'sans')+
  font("axis.text", size = 6,family = 'sans')+
  font("axis.title", size = 6,family = 'sans')+
  font("legend.title", size = 6,family = 'sans')+
  font("legend.text", size = 6,family = 'sans')+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))+  ####margin between top and bottom main and xy axis
  font("title", size = 6,family = 'sans')+
  font("subtitle", size = 6,family = 'sans')+
  rremove('xylab')+
  font('caption',size = 6,family = 'sans')+
  theme(strip.background =element_rect(fill = "transparent",  #####facet background
                                       colour = NA),
        panel.spacing = unit(0.1, "lines"))+        ####the distance between facet figures
  ggsave(filename = paste0('./','immune-score','.pdf'),useDingbats=FALSE,
         width = 170,height = 150,
         units = 'mm',dpi =600
  )

ICB.markers=read_xlsx(path =  "e:/immunogenemoci landsapes in subsites of HNSCC/immune checkpoint biomarkers.xlsx")


ICB.exp=mRNA[,ICB.markers$ICB.makers]%>%rownames_to_column(var = "ParticipantBarcode")



################# ICB heatmap and comparsions#######

meta<-clin%>%
  select(ParticipantBarcode,Subsite.1,Subsite.2,HPV_status)%>%merge(.,ICB.exp)


meta1=meta%>%filter(Subsite.2=="Non-OP")

colnames(meta1)[3]='Subsite'

meta1=meta1[,-2]

colnames(meta)[2]='Subsite'

meta=meta[,-3]

meta=rbind(meta,meta1)

meta$Subsite=factor(meta$Subsite,levels = c('Oropharynx',
                                            'Non-OP',
                                            'Hypopharynx',
                                            'Oral.cavity',
                                            'Larynx',
                                            'Cervical'))

meta=unite(meta,"HPV_status_Subsite",HPV_status,Subsite,remove = F)
#####################################ICB median compare and heatmap####################

t<-meta%>%split(.$HPV_status_Subsite)


tmp<-list()

for (i in names(t)) {tmp[[i]]=t[[i]]%>%select(5:36)}


Mean.dat.i<-sapply(tmp, function(x)colMeans(x,na.rm = T))


Mean.dat.i<-as.data.frame(Mean.dat.i)

mean.HPV<-t(Mean.dat.i)

mean.HPV<-as.data.frame(mean.HPV)
mean.HPV$Subsite<-rownames(mean.HPV)
mean.HPV<-mean.HPV%>%melt('Subsite',value.name = 'Proportion')
colnames(mean.HPV)[2]<-'Immune.cell'
table(mean.HPV$Subsite)

mean.HPV$HPV_status=stringr::str_split(mean.HPV$Subsite,"_",simplify = T)[,1]

mean.HPV$Subsite=stringr::str_split(mean.HPV$Subsite,"_",simplify = T)[,2]

mean.HPV$Subsite=as.factor(mean.HPV$Subsite)


mean.HPV$HPV_status=as.factor(mean.HPV$HPV_status)

mean.HPV$Subsite=factor(mean.HPV$Subsite,levels = c('Oropharynx',
                                                      'Non-OP',
                                                        'Hypopharynx',
                                                        'Oral.cavity',
                                                        'Larynx',
                                                        'Cervical'))


########


dl1=meta%>%select('Subsite','HPV_status',
                  ICB.markers$ICB.makers)%>%
  reshape2::melt(.,id.vars =c('Subsite',
                              'HPV_status'),
                 value.name = 'ICB.exp')%>%split(.,.$Subsite)
p=list()

foreach(i=1:6)%do%{
  
  a=compare_means(ICB.exp~HPV_status,data = dl1[[i]],method = 'wilcox.test',
                  group.by = 'variable',p.adjust.method = 'fdr')
  p[[i]]=a%>%select('variable','p.signif')
}  

names(p)=names(dl1)

p.mat=bind_rows(p, .id = "Subsite")

p.HPV=dcast(p.mat,Subsite~variable,value.var = 'p.signif')

rownames(p.HPV)=p.HPV$Subsite

p.HPV=p.HPV%>%select(-Subsite)

########################fold change##############

c2=list()

tmp2=mean.HPV%>%
  dcast(HPV_status+Subsite~Immune.cell,value.var = 'Proportion')%>%
  split(.,.$Subsite)

c2=lapply(names(tmp2), function(i){
  a=tmp2[[i]]%>%filter(HPV_status=='positive')%>%select(3:ncol(tmp2[[i]]))
  
  b=tmp2[[i]]%>%filter(HPV_status=='negative')%>%select(3:ncol(tmp2[[i]]))
  
  log2(a)-log2(b)})

names(c2)=names(tmp2)

tmp2=bind_rows(c2, .id = "Subsite")


rownames(tmp2)=tmp2$Subsite

tmp2=tmp2%>%select(-Subsite)

p.HPV=p.HPV[rownames(tmp2),]

p.HPV_1=p.HPV
p.HPV_1[p.HPV_1=="ns"]=0

library(circlize)
col_fun = colorRamp2(c(-5, 0,5), c("green", "white", "red"))
col_fun(seq(-5, 5))
cols=c("#00ba55",
       "#8a59ff",
       "#f5ff60",
       "#240084",
       "#092f00",
       "#9c0068",
       "#004362",
       "#ff9b96")
tmp2=as.matrix(tmp2)
ha = HeatmapAnnotation(ICB=colnames(tmp2),
                       annotation_name_gp = gpar(fontsize=6),
                       annotation_legend_param=list(
                         ICB =list(labels_gp = gpar(fontsize = 6),
                                   title_gp= gpar(fontsize = 6))
                       ),gp = gpar(col = "black")
)




pdf(file = "./heatmap.pdf",width = 12,height = 10)

Heatmap(tmp2,cluster_columns = F,cluster_rows = F,
        #col = col_fun,
        rect_gp = gpar(col = "white", lwd = 2),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", p.HPV[i, j]), x, y, gp = gpar(fontsize = 6))},
        
        top_annotation =ha,
        heatmap_width = unit(17, "cm"),
        heatmap_height = unit(6, "cm"),
        row_names_gp = grid::gpar(fontsize = 6),
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 6),
                                    title_gp= gpar(fontsize = 6)))

dev.off()
  

pdf(file = "./heatmap_1.pdf",width = 12,height = 10)

Heatmap(tmp2,cluster_columns = F,cluster_rows = F,
        #col = col_fun,
        rect_gp = gpar(col = "white", lwd = 2),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", p.HPV_1[i, j]), x, y, gp = gpar(fontsize = 6))},
        
        top_annotation =ha,
        heatmap_width = unit(17, "cm"),
        heatmap_height = unit(6, "cm"),
        row_names_gp = grid::gpar(fontsize = 6),
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 6),
                                    title_gp= gpar(fontsize = 6)))

dev.off()


