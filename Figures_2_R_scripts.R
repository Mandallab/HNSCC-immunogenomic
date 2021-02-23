rm(list = ls())

############set working file folder############
setwd(dir = "E:/immunogenemoci landsapes in subsites of HNSCC/manuscript/1.26/Figures_from_R/Figure 2/")
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
#options(scipen=6)

########Figure 2 KM curves######

load(file = 'E:/immunogenemoci landsapes in subsites of HNSCC/immunity score and genes/clin.data.HPV.HNSC.CESC.RData')

#CDR=readxl::read_xlsx(path = "F:/TCGA/TCGA-CDR-SupplementalTableS1.xlsx",sheet = 1)
#CDR=CDR%>%select(ParticipantBarcode=bcr_patient_barcode,OS,OS.time)

#clin.hpv_1=clin.hpv%>%select(-OS,-OS.time)


#df=merge(CDR,clin.hpv_1)


#df$OS.time=df$OS.time/30



#writexl::write_xlsx(df,path = "../../1.26/os_5_y.xlsx")

########univariate cox analysis#####

sur.mu.1=split(clin.hpv,clin.hpv$Subsite.2)


survfit(Surv(OS.time, 
                     OS==1)~HPV_status,data = sur.mu.1$Cervical)

summary(survfit(Surv(OS.time, 
                     OS==1)~HPV_status,sur.mu.1$Cervical),times = c(24,36)) #in months

summary(survfit(Surv(OS.time, 
                     OS==1)~HPV_status,sur.mu.1$Cervical),times = c(36)) #in months


summary(survfit(Surv(OS.time, 
                     OS==1)~HPV_status,sur.mu.1$OP),times = c(24,36)) #in months


summary(survfit(Surv(OS.time, 
                     OS==1)~HPV_status,sur.mu.1$`Non-OP`),times = c(24,36)) #in months






tmp=sur.mu.1[['Cervical']]%>%filter(HPV_status=='positive')#%>%
summary(survfit(Surv(OS.time, 
                     OS==1)~HPV_status,tmp),times = 60) #in months


tmp=sur.mu.1[['Cervical']]%>%filter(HPV_status=='negative')#%>%
summary(survfit(Surv(OS.time, 
                     OS==1)~HPV_status,tmp),times = 59) #in months


surv_summary(survfit(Surv(OS.time, 
                     OS==1)~HPV_status,tmp))





######follow-up time########

clin$OS.time[clin$OS.time==0]=NA


clin%>%select(ParticipantBarcode,Subsite.2,HPV_status,OS,OS.time,)%>%filter(Subsite.2=="OP")%>%summary.matrix()

clin%>%select(ParticipantBarcode,Subsite.2,HPV_status,OS,OS.time,)%>%filter(Subsite.2=="Non-OP")%>%summary.matrix()

clin%>%select(ParticipantBarcode,Subsite.2,HPV_status,OS,OS.time,)%>%filter(Subsite.2=="Cervical")%>%summary.matrix()



covariates=c("age","gender","Race",
             "Smoking",
             "T.stage",
             "N.stage",
             "M.stage",
             "Clinical.stage",
             "Histologic.grade",
             "Histologic.type",
             "HPV_status")

univ_formulas <- sapply(covariates,
                        function(x) {
                          as.formula(
                            paste('Surv(OS.time, 
                                  OS==1)~', x))
                        })

r=lapply(1:(length(sur.mu.1)-1), function(i){ 
  
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sur.mu.1[[i]])})
  
  # Extract data 
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$coefficients[,5], digits=3)
                           HR <-signif(x$coefficients[,2], digits=3);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                           HR <-paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-data.frame(HR, p.value)
                           rownames(res)<-rownames(x$coefficients)
                           return(res)
                         })
  res <- bind_rows(univ_results,.id = "column_label")
  
})

r[[3]]={c("age","Race",
          "Smoking",
          "T.stage",
          "N.stage",
          "M.stage",
          "Clinical.stage",
          "HPV_status")%>%
    sapply(.,function(x) {
      as.formula(
        paste('Surv(OS.time, 
                                  OS==1)~', x))
    })%>%lapply(., function(x){coxph(x, data = sur.mu.1[[3]])})%>% 
    lapply(.,
           function(x){ 
             x <- summary(x)
             p.value<-signif(x$coefficients[,5], digits=3)
             HR <-signif(x$coefficients[,2], digits=3);#exp(beta)
             HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
             HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
             HR <-paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
             res<-data.frame(HR, p.value)
             rownames(res)<-rownames(x$coefficients)
             return(res)
           })%>%bind_rows(.,.id = "column_label")
  
}

names(r)=names(sur.mu.1)



#############adjusted TNM age gender######OS###
######
{
  fit.1 <- coxph(Surv(OS.time, OS==1) ~ HPV_status+age+Smoking+gender+pathologic_T+pathologic_N+pathologic_M,
                 data=sur.mu.1$OP)
  
  fit.1%>%gtsummary::tbl_regression()
  
  
  p.1=paste0('p = ',round(summary(fit.1)[["coefficients"]][1,5],3))
  
  HR.1=paste0('HR = ',round(summary(fit.1)[["conf.int"]][1,1],2),
              ' (',
              round(summary(fit.1)[["conf.int"]][1,3],2),
              '-',
              round(summary(fit.1)[["conf.int"]][1,4],2),
              ')')
  ####
  
  
  fit.2 <- coxph(Surv(OS.time, OS==1) ~ HPV_status+age+Smoking+gender+pathologic_T+pathologic_N+pathologic_M,
                 data=sur.mu.1$`Non-OP`)
  
  fit.2%>%gtsummary::tbl_regression()
  
  
  p.2=paste0('p = ',round(summary(fit.2)[["coefficients"]][1,5],2))
  
  HR.2=paste0('HR = ',round(summary(fit.2)[["conf.int"]][1,1],2),
              ' (',
              round(summary(fit.2)[["conf.int"]][1,3],2),
              '-',
              round(summary(fit.2)[["conf.int"]][1,4],2),
              ')')
  
  
  fit.3 <- coxph(Surv(OS.time, OS==1) ~ HPV_status+age+Smoking+pathologic_T+pathologic_N+pathologic_M,
                 data=sur.mu.1$Cervical)
  
  fit.3%>%gtsummary::tbl_regression()
  
  p.3=paste0('p = ',round(summary(fit.3)[["coefficients"]][1,5],2))
  
  HR.3=paste0('HR = ',round(summary(fit.3)[["conf.int"]][1,1],2),
              ' (',
              round(summary(fit.3)[["conf.int"]][1,3],2),
              '-',
              round(summary(fit.3)[["conf.int"]][1,4],2),
              ')')
  
}
######################KM sur curve with risk table###################
xx=list()
foreach(i=1:3)%do%{
  p=ggsurvplot(surv_fit(Surv(OS.time, OS==1)~HPV_status,data = sur.mu.1[[i]]),
               risk.table = T, 
               linetype = c(1,2),
               risk.table.pos='out',
               #censor.size=2, size = 0.5,
               risk.table.fontsize=1.8,
               risk.table.y.text=F,
               legend.title='HPV_status',
               legend.labs=c('Negative','Positive'),
               pval.size=2.12,
               break.time.by=12,
               pval = T,
               palette='jco',
               tables.col = "strata",
               xlab = "Time (months)",
               title = paste0(LETTERS[i],'.Overall Survival in TCGA - ',names(sur.mu.1)[i],'(n=',nrow(sur.mu.1[[i]]),')'))
 # p$plot=p$plot+annotate(
 #   "text",
 #   x = 40, y = 0.25,
 #   vjust = 1, hjust = 1,
 #   label = paste0('Adjusted by age, gender and TNM stage','\n',
 #                  get(paste0('HR.',i)),'\n',get(paste0('p.',i))), 
 #   size=2.12
 # )
  p=p+theme_survminer(base_size = 6,
                      font.main = c(6, "bold", "darkblue"),
                      font.submain = c(6, "bold.italic", "purple"),
                      font.caption = c(6, "plain", "orange"),
                      font.x = c(6, "bold.italic", "black"),
                      font.y = c(6, "bold.italic", "darkred"),
                      font.tickslab = c(6, "plain", "darkgreen"),
                      font.legend = c(6,"bold",'black'))
  ggsave(print(p,newpage = F),filename = paste0('./KM-risktable-TCGA-OS- ',names(sur.mu.1)[i],'.pdf'),width = 90,height = 90,
         units = 'mm',dpi = 600)
  xx[[i]]=p
}  
  
  f2=arrange_ggsurvplots(xx, print = TRUE,
                         ncol = 3, nrow = 1, risk.table.height = 0.2)
  ggsave(f2,filename = './Fig2.Survival-OS.pdf',width = 200,height = 90,
         units = 'mm',dpi = 600)
  

  #########################matched cases for surv#########################
  
  
  library(MatchIt)
  
  dl<-clin
  dl[is.na(dl)]<-'unknown'
  dl<-na.omit(dl)
  
  dl$HPV_status.logic<-as.logical(dl$HPV_status=='positive')
  
  dl.1<-split(dl,dl$Subsite.2)
  m=list()
  set.seed(1234)
  for (x in names(dl.1)[1:2]) {
    
    match = matchit (HPV_status.logic ~HPV_status+age+gender+Smoking+pathologic_T+pathologic_N+pathologic_M,#age+Race+gender+pathologic_stage,
                     data =dl.1[[x]], method="nearest", ratio =1)
    df.match <- match.data(match)[1:ncol(dl.1[[x]])]
    
    m[[x]]<-df.match%>%select(ParticipantBarcode)
  }
  
  
  match = matchit (HPV_status.logic ~HPV_status+age+Smoking+pathologic_T+pathologic_N+pathologic_M,
                   data =dl.1[[3]], method="nearest", ratio =1)
  df.match <- match.data(match)[1:ncol(dl.1[[3]])]
  
  m[['Cervical']]<-df.match%>%select(ParticipantBarcode)
  
  m1<-bind_rows(m, .id = 'Subsite')
  
  #####################
  clin.matched<-subset(clin,clin$ParticipantBarcode%in%m1$ParticipantBarcode)

  clin.matched$Subsite.2=factor(clin.matched$Subsite.2,levels = c('OP','Non-OP','Cervical'))
  
  clin.matched$HPV_status[is.na(clin.matched$HPV_status)]<-'negative'
  
  clin.matched$HPV_status<-factor(clin.matched$HPV_status,levels = c('negative','positive'))

  #######################surivival #################
  
  
  ggsurvplot_facet(surv_fit(Surv(OS.time, OS==1)~HPV_status,data =  clin.matched),data =  clin.matched,facet.by = "Subsite.2",
                   palette = "jco", pval = TRUE, 
                   #surv.median.line = "hv",
                   break.time.by=12,
                   ncol = 3,censor.size=3,
  )+font("xy", size = 6)+font("axis.text", size = 6)+
    font("axis.title", size = 6)+
    font("legend.title", size = 6)+
    font("legend.text", size = 6)+
    ggsave(filename = './matched surv/matched_OS.pdf',width = 210,height = 210,
           units = 'mm',dpi = 300,
    )
  
  dl.3=split(clin.matched,clin.matched$Subsite.2)
  
  p1=ggsurvplot(surv_fit(Surv(OS.time, OS==1)~HPV_status,data = dl.3$OP),
               risk.table = 'nrisk_cumevents', 
               risk.table.fontsize=2,
               risk.table.y.text=F,
               pval.size=2,
               pval = T,
               break.time.by=12,
               palette='jco',
               tables.col = "strata",
               xlab = "time (months)",
               title = "A. Overall Survival in TCGA - paired OPSCC (n=50)")
  p1=p1+theme_survminer(base_size = 6,
                      
                      font.main = c(6, "bold", "darkblue"),
                      font.submain = c(6, "bold.italic", "purple"),
                      font.caption = c(6, "plain", "orange"),
                      font.x = c(6, "bold.italic", "black"),
                      font.y = c(6, "bold.italic", "darkred"),
                      font.tickslab = c(6, "plain", "darkgreen"),
                      font.legend = c(6,"bold",'black')
  )
  
  p1
  ggsave(print(p1,newpage = F),filename = './matched surv/OS-paired-OPSCC.pdf',width = 80,height = 90,
         units = 'mm',dpi = 300,
  )
  
  p2=ggsurvplot(surv_fit(Surv(OS.time, OS==1)~HPV_status,data = dl.3$`Non-OP`),
               risk.table = 'nrisk_cumevents', 
               risk.table.fontsize=2,
               break.time.by=12,
               risk.table.y.text=F,
               pval.size=2,
               pval = T,
               palette='jco',
               tables.col = "strata",
               xlab = "time (months)",
               title = "B. Overall Survival in TCGA - paired Non-OPSCC (n=54)")
  p2=p2+theme_survminer(base_size = 6,
                      
                      font.main = c(6, "bold", "darkblue"),
                      font.submain = c(6, "bold.italic", "purple"),
                      font.caption = c(6, "plain", "orange"),
                      font.x = c(6, "bold.italic", "black"),
                      font.y = c(6, "bold.italic", "darkred"),
                      font.tickslab = c(6, "plain", "darkgreen"),
                      font.legend = c(6,"bold",'black')
  )
  
  p2
  ggsave(print(p2,newpage = F),filename = './matched surv/OS-paired-Non-OPSCC.pdf',width = 80,height = 90,
         units = 'mm',dpi = 300,
  )
  
  p3=ggsurvplot(surv_fit(Surv(OS.time, OS==1)~HPV_status,data = dl.3$Cervical),
               risk.table = 'nrisk_cumevents', 
               risk.table.fontsize=2,
               risk.table.y.text=F,
               break.time.by=12,
               pval.size=2,
               pval = T,
               palette='jco',
               tables.col = "strata",
               xlab = "time (months)",
               title = "C. Overall Survival in TCGA - paired CESC (n=20)")
  p3=p3+theme_survminer(base_size = 6,
                      
                      font.main = c(6, "bold", "darkblue"),
                      font.submain = c(6, "bold.italic", "purple"),
                      font.caption = c(6, "plain", "orange"),
                      font.x = c(6, "bold.italic", "black"),
                      font.y = c(6, "bold.italic", "darkred"),
                      font.tickslab = c(6, "plain", "darkgreen"),
                      font.legend = c(6,"bold",'black')
  )
  
  p3
  ggsave(print(p3,newpage = F),filename = './matched surv/OS-paired-CESC.pdf',width = 80,height = 90,
         units = 'mm',dpi = 300,
  )

yy=list()
  

yy[[1]]=p1

yy[[2]]=p2

yy[[3]]=p3

ff2=arrange_ggsurvplots(yy,ncol = 3)

ggsave(ff2,filename = './matched surv/Survival-OS.pdf',width = 200,height = 90,
units = 'mm',dpi = 600)   

fit=surv_fit(Surv(OS.time, OS==1) ~ HPV_status,
      data= dl.3$OP)

surv_pvalue(fit,dl.3$OP)

survdiff(Surv(OS.time, OS==1) ~ HPV_status,
         data= dl.3$OP)

############

dl.3$OP%>%select(age,Smoking,gender,pathologic_T,
                 pathologic_N,pathologic_M,HPV_status,
                 Subsite.2)%>%
  tbl_summary(by=HPV_status)%>%add_p()

dl.3$`Non-OP`%>%select(age,Smoking,gender,pathologic_T,
                 pathologic_N,pathologic_M,HPV_status,
                 Subsite.2)%>%
  tbl_summary(by=HPV_status)%>%add_p()


dl.3$Cervical%>%select(age,Smoking,pathologic_T,
                 pathologic_N,pathologic_M,HPV_status,
                 Subsite.2)%>%
  tbl_summary(by=HPV_status)%>%add_p()


clin.matched%>%
  table1::table1(~age+Smoking+pathologic_T+pathologic_N+pathologic_M|HPV_status+Subsite.2)
