pkg_control(c('ComplexHeatmap','stringr',"maftools","tidyverse","ggplot2",'psych','rcompanion','mosaicData','ggcorrplot','ggplot2','msigdbr','annotate','org.Hs.eg.db',
              'factoextra','survminer','readxl','TiMEx','glmnet','vcdExtra','glmnet','ggmosaic','survival','circlize','tm','GSgalgoR','lsa','arkhe','wordcloud'))

###################################################################################################
####################################  Genetic profile summary  ####################################
###################################################################################################
top_genes=c('KRAS','TP53','CDKN2A','SMAD4','MSH3','ARID1A','SPTA1','BRCA2','KMT2C','RNF43')
genes <- c(top_genes)
FMG_mat <- pdac_mat_filled[,paste0(c(both_sample,FFX_sample,AG_sample,gray_sample),'D')]; FMG_mat <- FMG_mat[rowSums(FMG_mat == "") != ncol(FMG_mat),]
FMGs <- genes
FMG_mat[FMG_mat == ""] <- "background"



'AGE'=na.omit(meta$AGE[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))])
'TUMORSIZE'=na.omit(meta$TUMORSIZE[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))])
'MSI'=na.omit(meta$MSI[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))])
'TMB'=na.omit(meta$TMB[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))])

top_anno <- HeatmapAnnotation(
  'SEX'=factor(ifelse(meta$SEX[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))] == 1, "M", "F"),levels=c('F','M')),
  'LOCATION'=factor(ifelse(meta$LOCATION[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))] == 0, "HEAD", ifelse(meta$LOCATION[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))] == 1, "BODY", ifelse(meta$LOCATION[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))] == 2, "TAIL", ifelse(meta$LOCATION[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))] == 3, "MIXED", NA)))),levels=c('HEAD','BODY','TAIL','MIXED','NA')),
  'AGE'=meta$AGE[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))],
  'TUMORSIZE'=meta$TUMORSIZE[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))],
  'MSI'=meta$MSI[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))],
  'TMB'=meta$TMB[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))],
  'FFX BR'=factor(ifelse(meta$FFX_BR[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))] == 0, "CR/PR", ifelse(meta$FFX_BR[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))] == 1, "SD", ifelse(meta$FFX_BR[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))] == 2, "PD", NA))),levels=c('CR/PR','SD','PD','NA')),
  'AG BR'=factor(ifelse(meta$AG_BR[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))] == 0, "CR/PR", ifelse(meta$AG_BR[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))] == 1, "SD", ifelse(meta$AG_BR[match(colnames(FMG_mat),paste0(meta$SAMPLE,'D'))] == 2, "PD", NA))),levels=c('CR/PR','SD','PD','NA')),
  "Response" = factor(c(rep("BOTH", length(both_sample)), rep("FFX ONLY", length(FFX_sample)), rep("AG ONLY", length(AG_sample)), rep("NONE / GRAY ZONE", length(gray_sample))),levels = c("BOTH","FFX ONLY","AG ONLY","NONE / GRAY ZONE")),
  
  col = list("SEX" = c("M" ="#EAC7C7", "F" ="#A0C3D2"),
             "LOCATION" = c("HEAD" ="#C4DFDF", "BODY" ="#D2E9E9", "TAIL" ="#E3F4F4", "MIXED" ="#F8F6F4", "NA" ="#F8F6F4"),
             "AGE" = colorRamp2(c(min(AGE), max(AGE)), c("white", "#8D8DAA")),
             "TUMORSIZE" = colorRamp2(c(min(TUMORSIZE), max(TUMORSIZE)), c("white", "#8D8DAA")),
             "MSI" = colorRamp2(c(min(MSI), max(MSI)), c("white", "#8D8DAA")),
             "TMB" = colorRamp2(c(min(TMB), max(TMB)), c("white", "#8D8DAA")),
             "FFX BR" = c("CR/PR" ="#85A389", "SD" ="#A2CDB0", "PD" ="#F1C27B", "NA" ="#F8F6F4"),
             "AG BR" = c("CR/PR" ="#85A389", "SD" ="#A2CDB0", "PD" ="#F1C27B", "NA" ="#F8F6F4"),
             "Response" = c("BOTH" ="#EB4747",  "FFX ONLY" ="#FF8B8B", "AG ONLY" ="#FFDEDE","NONE / GRAY ZONE" ="#ABC9FF")),
  show_annotation_name = TRUE,
  simple_anno_size = unit(4, "mm")
)



alter_fun$DUP <- function(x, y, w, h) {grid.rect(x, y, w * 0.8, h * 0.75, gp = gpar(lwd = 2, col = '#EB6440', fill = NA))}
oncoPrint(FMG_mat[FMGs, paste0(c(both_sample,FFX_sample,AG_sample,gray_sample), 'D')],
          row_order = FMGs,
          column_split = factor(c(rep("BOTH", length(both_sample)), rep("FFX ONLY", length(FFX_sample)), rep("AG ONLY", length(AG_sample)), rep("NONE / GRAY ZONE", length(gray_sample))),levels = c("BOTH","FFX ONLY","AG ONLY","NONE / GRAY ZONE")),
          alter_fun = alter_fun,
          top_annotation = top_anno,
          remove_empty_columns = TRUE,
          remove_empty_rows = FALSE,
          show_row_names = TRUE,
          show_column_names = FALSE,
          right_annotation = NULL,
          pct_side = "left",
          show_pct = FALSE,
          row_names_side = "right",
          alter_fun_is_vectorized = FALSE,
          show_heatmap_legend = TRUE)


###################################################################################################
################################  Clinical pathogenic features  ###################################
###################################################################################################
meta_test=METADATA[METADATA$FFX %in% c(1,2,3),]
meta_test$SEX=factor(meta_test$SEX,levels = c(0,1))
meta_test$MSI=as.numeric(meta_test$MSI)
meta_test$MSI <- ifelse(is.na(meta_test$MSI),  NA, ifelse(meta_test$MSI < 3.5,  "MSS", ifelse(meta_test$MSI < 10,  "MSI-L", meta_test$MSI)))
meta_test$MSI <- factor(meta_test$MSI, levels = c("MSS", "MSI-L"))
meta_test$TMB  <- ifelse(meta_test$TMB >=mean(na.omit(meta$TMB)), "TMB high", "TMB low")
meta_test$TMB=as.factor(meta_test$TMB)
meta_test$LOCATION=factor(meta_test$LOCATION,levels = c(0,1,2))
meta_test$STAGE <- factor(ifelse(meta_test$STAGE == 1, "M0", "M1"))
res_tmp<-survival::coxph(survival::Surv(FFX_PFS)~SEX+AGE+MSI+LOCATION+TMB+TUMORSIZE+STAGE,data=meta_test) 
a=ggforest(res_tmp,data=meta_test,main = "FFX");a



meta_test=METADATA[METADATA$AG %in% c(1,2,3),]
meta_test$SEX=factor(meta_test$SEX,levels = c(0,1))
meta_test$MSI=as.numeric(meta_test$MSI)
meta_test$MSI <- ifelse(is.na(meta_test$MSI),  NA, ifelse(meta_test$MSI < 3.5,  "MSS", ifelse(meta_test$MSI < 10,  "MSI-L", meta_test$MSI)))
meta_test$MSI <- factor(meta_test$MSI, levels = c("MSS", "MSI-L"))
meta_test$TMB  <- ifelse(meta_test$TMB >=mean(na.omit(meta$TMB)), "TMB high", "TMB low")
meta_test$TMB=as.factor(meta_test$TMB)
meta_test$LOCATION=factor(meta_test$LOCATION,levels = c(0,1,2))
meta_test$STAGE <- factor(ifelse(meta_test$STAGE == 1, "M0", "M1"))
res_tmp<-survival::coxph(survival::Surv(AG_PFS)~SEX+AGE+MSI+LOCATION+TMB+TUMORSIZE+STAGE,data=meta_test)
b=ggforest(res_tmp, data = meta_test,main = "AG");b
a|b
