pkg_control(c('ComplexHeatmap','stringr',"maftools","tidyverse","ggplot2",'psych','rcompanion','mosaicData','ggcorrplot','ggplot2','msigdbr','annotate','org.Hs.eg.db',
              'factoextra','survminer','readxl','TiMEx','glmnet','vcdExtra','glmnet','ggmosaic','survival','circlize','tm','GSgalgoR','lsa','arkhe','wordcloud'))

#########################################################################################
###########################  DDR+TGFbeta pathway mutation  ##############################
#########################################################################################
meta <- meta_YUHS
PDAC_lib_df <- PDAC_lib_df_YUHS
PDAC_MAFs <- PDAC_MAFs_YUHS 
PDAC_VAF_mat <- PDAC_VAF_mat_norm_YUHS 
ddr_gene1 <- c('MLH1','MSH2','MSH6','PMS1','PMS2','EPCAM')
ddr_gene2 <- c('ERCC2','ERCC3','ERCC4','ERCC5')
ddr_gene3 <- c('BRCA1','MRE11A','NBN','RAD50','RAD51','RAD51B','RAD51D','RAD52','RAD54L')
ddr_gene4 <- c('BRCA2','BRIP1','FANCA','FANCC','PALB2','RAD51C','BLM')
ddr_gene5 <- c('ATM','ATR','CHEK1','CHEK2','MDC1')
ddr_gene6 <- c('POLE','MUTYH','PARP1','RECQL4')
ddr_gene <- c(ddr_gene1,ddr_gene2,ddr_gene3,ddr_gene4,ddr_gene5)
tso_gene <- read.table("/home/eunu0605/Project/PDAC/Library/GeneList/TSO_genelist.txt"); tso_gene <- tso_gene$V1
TGF_gene=c('SMAD4','SMAD3','TGFBR1','TGFBR2','INHBA','ACVR1B','CBL','MAPK3','MEN1','MYC','NCOR1','PARP1','XPO1')

let=unique(c(TGF_gene,ddr_gene))
let <- let[let %in% colnames(meta)]
let=let[let %in% colnames(meta_YUHS)][colSums(meta_YUHS[meta_YUHS$FFX %in% c(1,2,3,4) ,let[let %in% colnames(meta_YUHS)]])>=2]; length(let)
let=c(let[!(let %in% c( "SMAD4") )],'SMAD4type')
SMAD4_up_sample=c() ; SMAD4_down_sample=c()
for(i in let){
  if(i=='SMAD4type'){
    SMAD4_mut_sample=meta$SAMPLE[meta$SMAD4type==1]
    down_sample=SMAD4_mut_sample[PDAC_VAF_mat['SMAD4',SMAD4_mut_sample]<0.2]
    up_sample=SMAD4_mut_sample[PDAC_VAF_mat['SMAD4',SMAD4_mut_sample]>=0.2]
  }else {
    SMAD4_mut_sample=colnames(PDAC_VAF_mat)[PDAC_VAF_mat[i,]!=0]
    down_sample=SMAD4_mut_sample[PDAC_VAF_mat[i,SMAD4_mut_sample]<0.2]
    up_sample=SMAD4_mut_sample[PDAC_VAF_mat[i,SMAD4_mut_sample]>=0.2 ]
  }
  print(i)
  SMAD4_up_sample=c(SMAD4_up_sample,up_sample)
  SMAD4_down_sample=c(SMAD4_down_sample,down_sample)
}
check_meta <- meta[meta$FFX %in% c(1,2,3),]
GOI=let ; GOI <- GOI[GOI %in% colnames(check_meta)] 
check_meta$group <- rowSums(check_meta[,GOI[colSums(check_meta[,GOI])>0]]);check_meta$group
check_meta$group[check_meta$group>0] <- 'MUT' ;check_meta$group[check_meta$group<=0]<- 'WT'
check_meta$group <- factor(check_meta$group)
fit<-survfit(Surv(FFX_PFS)~group, data=check_meta)
YUHS_res3=survdiff(Surv(FFX_PFS)~group, data=check_meta)$p
print(YUHS_res3)
ggsurv <- ggsurvplot(fit, data=check_meta, #pval = TRUE,
                     palette = c("#F45050", "#7C73C0"),
                     xlim = c(0,500), xlab = "Time(Days)",
                     break.time.by = 100,
                     legend.title=''
);ggsurv
ggsurv$plot <- ggsurv$plot + theme(
  legend.text = element_text(size = 15, face = "bold"),
  legend.title = element_text(size = 15, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
) + labs(title = "TGF+DDR geneset (p=0.003)"); ggsurv

###################################################################################################################
###########################  DDR+TGFbeta pathway mutation (first-line FFX treatment) ##############################
###################################################################################################################
let=unique(c(TGF_gene,ddr_gene))
let <- let[let %in% colnames(meta)]
let=let[let %in% colnames(meta_YUHS)][colSums(meta_YUHS[meta_YUHS$FFX %in% c(1,2,3,4) ,let[let %in% colnames(meta_YUHS)]])>=2]; length(let)
let=c(let[!(let %in% c( "SMAD4") )],'SMAD4type')
check_meta <- meta[meta$FFX %in% c(1,2,3) & !is.na(meta$first_FFX),]
GOI=let ; GOI <- GOI[GOI %in% colnames(check_meta)] 
check_meta$group <- rowSums(check_meta[,GOI[colSums(check_meta[,GOI])>0]]);check_meta$group
check_meta$group[check_meta$group>0] <- 'MUT' ;check_meta$group[check_meta$group<=0]<- 'WT'
check_meta$group <- factor(check_meta$group)
fit<-survfit(Surv(FFX_PFS)~group, data=check_meta)
YUHS_res3=survdiff(Surv(FFX_PFS)~group, data=check_meta)$p
print(YUHS_res3)
ggsurv <- ggsurvplot(fit, data=check_meta,
                     palette = c("#F45050", "#7C73C0"),
                     xlim = c(0,500), xlab = "Time(Days)",
                     break.time.by = 100,
                     legend.title=''
);ggsurv
ggsurv$plot <- ggsurv$plot + theme(
  legend.text = element_text(size = 15, face = "bold"),
  legend.title = element_text(size = 15, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
) + labs(title = "TGF+DDR geneset (p=0.02)"); ggsurv

#############################################################################################################
###########################  DDR+TGFbeta pathway mutation (FFX only treatment) ##############################
#############################################################################################################
let=unique(c(TGF_gene,ddr_gene))
let <- let[let %in% colnames(meta)]
let=let[let %in% colnames(meta_YUHS)][colSums(meta_YUHS[meta_YUHS$FFX %in% c(1,2,3,4) ,let[let %in% colnames(meta_YUHS)]])>=2]; length(let)
let=c(let[!(let %in% c( "SMAD4") )],'SMAD4type')
check_meta <- meta[meta$FFX %in% c(1,2,3) & meta$AG==4,]
GOI=let ; GOI <- GOI[GOI %in% colnames(check_meta)] 
check_meta$group <- rowSums(check_meta[,GOI[colSums(check_meta[,GOI])>0]]);check_meta$group
check_meta$group[check_meta$group>0] <- 'MUT' ;check_meta$group[check_meta$group<=0]<- 'WT'
check_meta$group <- factor(check_meta$group)
fit<-survfit(Surv(FFX_PFS)~group, data=check_meta)
YUHS_res3=survdiff(Surv(FFX_PFS)~group, data=check_meta)$p
print(YUHS_res3)
ggsurv <- ggsurvplot(fit, data=check_meta, 
                     palette = c("#F45050", "#7C73C0"),
                     xlim = c(0,500), xlab = "Time(Days)",
                     break.time.by = 100,
                     legend.title=''
);ggsurv
ggsurv$plot <- ggsurv$plot + theme(
  legend.text = element_text(size = 15, face = "bold"),
  legend.title = element_text(size = 15, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
) + labs(title = "TGF+DDR geneset (p=0.04)"); ggsurv

#########################################################################################################
###########################  DDR+TGFbeta pathway mutation (high clonality) ##############################
#########################################################################################################
check_meta <- meta[meta$FFX %in% c(1,2,3),] 
check_meta$group <- 0 ; check_meta$group[check_meta$SAMPLE %in% c(SMAD4_down_sample)] <- 0 
check_meta$group[check_meta$SAMPLE %in% SMAD4_up_sample] <-1
check_meta <- check_meta[check_meta$FFX_PFS!=0,]
check_meta$group <- factor(check_meta$group)
fit<-survfit(Surv(FFX_PFS)~group, data=check_meta, conf.lower='modified')
diff<-survdiff(Surv(FFX_PFS)~group, data=check_meta) ; diff$pvalue
ggsurv <- ggsurvplot(fit, data=check_meta, 
                     palette = c("#7C73C0", "#F45050"),
                     xlim = c(0,500), xlab = "Time(Days)",
                     break.time.by = 100,
                     legend.labs = c("MUT=NO", "MUT=YES"),
                     legend.title=''
);ggsurv
ggsurv$plot <- ggsurv$plot + theme(
  legend.text = element_text(size = 15, face = "bold"),
  legend.title = element_text(size = 15, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
) + labs(title = paste0("TGF+DDR High VAF (p=0.001)")); ggsurv