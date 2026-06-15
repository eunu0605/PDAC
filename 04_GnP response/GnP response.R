pkg_control(c('ComplexHeatmap','stringr',"maftools","tidyverse","ggplot2",'psych','rcompanion','mosaicData','ggcorrplot','ggplot2','msigdbr','annotate','org.Hs.eg.db',
              'factoextra','survminer','readxl','TiMEx','glmnet','vcdExtra','glmnet','ggmosaic','survival','circlize','tm','GSgalgoR','lsa','arkhe','wordcloud'))

#################################################################################
###########################  p53 pathway mutation  ##############################
#################################################################################
meta <- meta_YUHS
PDAC_lib_df <- PDAC_lib_df_YUHS
PDAC_MAFs <- PDAC_MAFs_YUHS 
PDAC_VAF_mat <- PDAC_VAF_mat_norm_YUHS 
p53_gene <-c('ATM','ATR','MDC1','RAD51D','MSH2','MLH1','FANCI','FANCD2','FANCC','BRCA1','TP53')
p53_gene=unique(c(p53_gene))
p53_gene <- p53_gene[p53_gene %in% colnames(meta)]
p53_gene=p53_gene[p53_gene %in% colnames(meta_YUHS)][colSums(meta_YUHS[meta_YUHS$AG %in% c(1,2,3,4) ,p53_gene[p53_gene %in% colnames(meta_YUHS)]])>=2]; length(p53_gene)
let=p53_gene
let=c(let[!(let %in% c( 'TP53') )],'TP53type')
SMAD4_up_sample=c() ; SMAD4_down_sample=c()
for(i in let){
  if(i=='TP53type'){
    SMAD4_mut_sample=meta$SAMPLE[meta$TP53type==1]
    down_sample=SMAD4_mut_sample[PDAC_VAF_mat['TP53',SMAD4_mut_sample]<0.2]
    up_sample=SMAD4_mut_sample
    up_sample=SMAD4_mut_sample[PDAC_VAF_mat['TP53',SMAD4_mut_sample]>=0.2 ]
  }else if(i=='CDKN2Atype'){
    SMAD4_mut_sample=meta$SAMPLE[meta$CDKN2Atype==1]
    down_sample=SMAD4_mut_sample[PDAC_VAF_mat['CDKN2A',SMAD4_mut_sample]<0.2]
    up_sample=SMAD4_mut_sample
    up_sample=SMAD4_mut_sample[PDAC_VAF_mat['CDKN2A',SMAD4_mut_sample]>=0.2 ]
  }else{
    SMAD4_mut_sample=colnames(PDAC_VAF_mat)[PDAC_VAF_mat[i,]!=0]
    down_sample=SMAD4_mut_sample[PDAC_VAF_mat[i,SMAD4_mut_sample]<0.2]
    up_sample=SMAD4_mut_sample
    up_sample=SMAD4_mut_sample[PDAC_VAF_mat[i,SMAD4_mut_sample]>=0.2 ]
  }
  SMAD4_up_sample=c(SMAD4_up_sample,up_sample)
  SMAD4_down_sample=c(SMAD4_down_sample,down_sample)
  print(i)
}
check_meta <- meta[meta$AG %in% c(1,2,3),] 
GOI=let; GOI <- GOI[GOI %in% colnames(meta)] ; GOI <- c(GOI)
check_meta$group <- rowSums(check_meta[,GOI[colSums(check_meta[,GOI])>0]]);check_meta$group 
check_meta$group[check_meta$group>0] <- 'MUT' ;check_meta$group[check_meta$group<=0]<- 'WT'

check_meta$group <- factor(check_meta$group)
fit<-survfit(Surv(AG_PFS)~group, data=check_meta)
survdiff(Surv(AG_PFS)~group, data=check_meta)
ggsurv <- ggsurvplot(fit, data=check_meta, 
                     palette = c("#A8D1D1", "#F6D776"),
                     xlim = c(0,500), xlab = "Time(Days)",
                     break.time.by = 100,
                     legend.title='',
                     legend.labs = c("MUT=YES", "MUT=NO")
);ggsurv
ggsurv$plot <- ggsurv$plot + theme(
  legend.text = element_text(size = 15, face = "bold"),
  legend.title = element_text(size = 15, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
) + labs(title = "p53 geneset (p=0.02)"); ggsurv

###########################################################################################################
###########################  p53 pathway mutation (first-line Gnp treatment) ##############################
###########################################################################################################
check_meta <- meta[meta$AG %in% c(1,2,3) &!is.na(meta$first_AG),] 
GOI=let; GOI <- GOI[GOI %in% colnames(meta)] ; GOI <- c(GOI)
check_meta$group <- rowSums(check_meta[,c(GOI[colSums(check_meta[,GOI])>0],GOI[colSums(check_meta[,GOI])>0])]);check_meta$group 
check_meta$group[check_meta$group>0] <- 'MUT' ;check_meta$group[check_meta$group<=0]<- 'WT'

check_meta$group <- factor(check_meta$group)
fit<-survfit(Surv(AG_PFS)~group, data=check_meta)
survdiff(Surv(AG_PFS)~group, data=check_meta)
ggsurv <- ggsurvplot(fit, data=check_meta, #pval = TRUE,
                     palette = c("#A8D1D1", "#F6D776"),
                     xlim = c(0,500), xlab = "Time(Days)",
                     break.time.by = 100,
                     legend.title='',
                     legend.labs = c("MUT=YES", "MUT=NO")
);ggsurv
ggsurv$plot <- ggsurv$plot + theme(
  legend.text = element_text(size = 15, face = "bold"),
  legend.title = element_text(size = 15, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
) + labs(title = "p53 geneset (p=0.02)"); ggsurv

#####################################################################################################
###########################  p53 pathway mutation (GnP only treatment) ##############################
#####################################################################################################
check_meta <- meta[meta$AG %in% c(1,2,3) &meta$FFX==4,] 
GOI=let; GOI <- GOI[GOI %in% colnames(meta)] ; GOI <- c(GOI)
check_meta$group <- rowSums(check_meta[,c(GOI[colSums(check_meta[,GOI])>0],GOI[colSums(check_meta[,GOI])>0])]);check_meta$group 
check_meta$group[check_meta$group>0] <- 'MUT' ;check_meta$group[check_meta$group<=0]<- 'WT'

check_meta$group <- factor(check_meta$group)
fit<-survfit(Surv(AG_PFS)~group, data=check_meta)
survdiff(Surv(AG_PFS)~group, data=check_meta)
ggsurv <- ggsurvplot(fit, data=check_meta, 
                     palette = c("#A8D1D1", "#F6D776"),
                     xlim = c(0,500), xlab = "Time(Days)",
                     break.time.by = 100,
                     legend.title='',
                     legend.labs = c("MUT=YES", "MUT=NO")
);ggsurv
ggsurv$plot <- ggsurv$plot + theme(
  legend.text = element_text(size = 15, face = "bold"),
  legend.title = element_text(size = 15, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
) + labs(title = "p53 geneset (p=0.0005)"); ggsurv

#################################################################################################
###########################  p53 pathway mutation (high clonality) ##############################
##################################################################################################
check_meta <- meta[meta$AG %in% c(1,2,3),] 
check_meta$group <- 0 ; check_meta$group[check_meta$SAMPLE %in% c(SMAD4_down_sample)] <- 0 ; 
check_meta$group[check_meta$SAMPLE %in% SMAD4_up_sample] <-1
check_meta$group <- factor(check_meta$group)
check_meta <- check_meta[check_meta$AG_PFS!=0,]
fit<-survfit(Surv(AG_PFS)~group, data=check_meta, conf.lower='modified')
diff<-survdiff(Surv(AG_PFS)~group, data=check_meta); print(diff)
ggsurv <- ggsurvplot(fit, data=check_meta, 
                     palette = c("#F6D776", "#A8D1D1"),
                     xlim = c(0,500), xlab = "Time(Days)",
                     break.time.by = 100,
                     legend.labs = c("MUT=NO", "MUT=YES"),
                     legend.title=''
);ggsurv
ggsurv$plot <- ggsurv$plot + theme(
  legend.text = element_text(size = 15, face = "bold"),
  legend.title = element_text(size = 15, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
) + labs(title = paste0("p53 High VAF (p=0.0009)")); ggsurv
