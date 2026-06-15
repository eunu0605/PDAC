pkg_control(c('ComplexHeatmap','stringr',"maftools","tidyverse","ggplot2",'psych','rcompanion','mosaicData','ggcorrplot','ggplot2','msigdbr','annotate','org.Hs.eg.db',
              'factoextra','survminer','readxl','TiMEx','glmnet','vcdExtra','glmnet','ggmosaic','survival','circlize','tm','GSgalgoR','lsa','arkhe','wordcloud'))

###################################################################################################
###########################  FOLFIRINOX response related major gene  ##############################
###################################################################################################
rFFX_mut_mat <- pdac_mut_mat_bin[paste0(rFFX,'D'),]; rFFX_mut_mat <- rFFX_mut_mat[,colSums(rFFX_mut_mat)!=0]
mrFFX_mut_mat <- pdac_mut_mat_bin[paste0(mrFFX,'D'),]; mrFFX_mut_mat <- mrFFX_mut_mat[,colSums(mrFFX_mut_mat)!=0]
nrFFX_mut_mat <- pdac_mut_mat_bin[paste0(nrFFX,'D'),]; nrFFX_mut_mat <- nrFFX_mut_mat[,colSums(nrFFX_mut_mat)!=0]
top_gene=colnames(pdac_mut_mat_bin)[colSums(pdac_mut_mat_bin)>(length(rFFX)+length(nrFFX))*0.25]

Onco_mat <- data.frame(matrix(NA,
                              length(top_gene),3,
                              dimnames = list(top_gene,
                                              c('rFFX','mrFFX','nrFFX'))))
for(group in c('rFFX','mrFFX','nrFFX')){
  mut_mat <- get(paste0(group,'_mut_mat'))
  for(idx in seq(length(top_gene))){
    gene <-top_gene[idx]
    if(gene %in% colnames(mut_mat)){
      value <- sum(mut_mat[,gene])
    } else{value <- 0}
    Onco_mat[idx,group] <- value
  }
};check.del <- rowSums(Onco_mat)!=0; Onco_mat <-  Onco_mat[check.del,];Onco_mat

fish_g_pval <- c() 
for(idx in seq(nrow(Onco_mat))){
  tmp_table <- matrix(c(Onco_mat[idx,1],length(rFFX)-Onco_mat[idx,1], Onco_mat[idx,3], length(nrFFX)-Onco_mat[idx,3] ) , nrow = 2)
  fisher <- fisher.test(tmp_table) 
  fish_g_pval <- c(fish_g_pval,fisher$p.value) 
} ; Onco_mat$pval <- fish_g_pval
Onco_mat[sapply(strsplit(rownames(Onco_mat), '\\.'), tail, n = 1) %in% 'SMAD4',]
Onco_mat[Onco_mat$pval<0.05 ,]

fish_adj <- p.adjust(Onco_mat$pval,method='fdr')
Onco_mat$adj <- fish_adj
Onco_mat

total_num=nrow(rFFX_mut_mat)+nrow(nrFFX_mut_mat)
folfirinox <- data.frame(MUT=c(rep('MUT',Onco_mat['SMAD4',1]),rep('WT',nrow(rFFX_mut_mat)-Onco_mat['SMAD4',1]),rep('MUT',Onco_mat['SMAD4',3]),rep('WT',nrow(nrFFX_mut_mat)-Onco_mat['SMAD4',3])),
                         RES=c(rep('RES',nrow(rFFX_mut_mat)),rep('NONRES',nrow(nrFFX_mut_mat))))
folfirinox %>%
  ggplot(  ) +
  geom_mosaic(aes(x = product(RES), fill = MUT)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  ylab("SMAD4 mutation frequency") +
  scale_fill_manual(values = c("#F45050", "grey90")) +  
  theme_bw() +
  theme(axis.title.y = element_text(size = 14, face = "bold"))


mosaicplot(~ RES+MUT, data = folfirinox, color = TRUE, cex = 1,.2)


###################################################################################################
##################  Gemcitabine nab paclitaxel response related major gene  #######################
###################################################################################################
rAG_mut_mat <- pdac_mut_mat_bin[paste0(rAG,'D'),]; rAG_mut_mat <- rAG_mut_mat[,colSums(rAG_mut_mat)!=0]
mrAG_mut_mat <- pdac_mut_mat_bin[paste0(mrAG,'D'),]; mrAG_mut_mat <- mrAG_mut_mat[,colSums(mrAG_mut_mat)!=0]
nrAG_mut_mat <- pdac_mut_mat_bin[paste0(nrAG,'D'),]; nrAG_mut_mat <- nrAG_mut_mat[,colSums(nrAG_mut_mat)!=0]
top_gene=colnames(pdac_mut_mat_bin)[colSums(pdac_mut_mat_bin)>(length(rAG)+length(nrAG))*0.25]

Onco_mat <- data.frame(matrix(NA,
                              length(top_gene),3,
                              dimnames = list(top_gene,
                                              c('rAG','mrAG','nrAG'))))
for(group in c('rAG','mrAG','nrAG')){
  mut_mat <- get(paste0(group,'_mut_mat'))
  for(idx in seq(length(top_gene))){
    gene <-top_gene[idx]
    if(gene %in% colnames(mut_mat)){
      value <- sum(mut_mat[,gene])
    } else{value <- 0}
    Onco_mat[idx,group] <- value
  }
};check.del <- rowSums(Onco_mat)!=0; Onco_mat <-  Onco_mat[check.del,];Onco_mat

fish_g_pval <- c() 
for(idx in seq(nrow(Onco_mat))){
  tmp_table <- matrix(c(Onco_mat[idx,1],length(rFFX)-Onco_mat[idx,1], Onco_mat[idx,3], length(nrFFX)-Onco_mat[idx,3] ) , nrow = 2)
  fisher <- fisher.test(tmp_table) 
  fish_g_pval <- c(fish_g_pval,fisher$p.value) 
} ; Onco_mat$pval <- fish_g_pval
Onco_mat[sapply(strsplit(rownames(Onco_mat), '\\.'), tail, n = 1) %in% 'SMAD4',]
Onco_mat[Onco_mat$pval<0.05 ,]

fish_adj <- p.adjust(Onco_mat$pval,method='fdr')
Onco_mat$adj <- fish_adj
Onco_mat

total_num=nrow(rAG_mut_mat)+nrow(nrAG_mut_mat)
gem <- data.frame(MUT=c(rep('MUT',Onco_mat['TP53',1]),rep('WT',nrow(rAG_mut_mat)-Onco_mat['TP53',1]),rep('MUT',Onco_mat['TP53',3]),rep('WT',nrow(nrAG_mut_mat)-Onco_mat['TP53',3])),
                  RES=c(rep('RES',nrow(rAG_mut_mat)),rep('NONRES',nrow(nrAG_mut_mat))))
gem %>%
  ggplot() +
  geom_mosaic(aes(x = product(RES), fill = MUT)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(values = c("#A8D1D1", "grey90")) +
  ylab("TP53 mutation frequency") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14, face = "bold"))

###################################################################################################
###############  FOLFIRINOX response related major gene : SMAD4 lollipopPlot  #####################
###################################################################################################
lollipopPlot(
  maf = PDAC_MAFs,
  gene = 'SMAD4',
  AACol= 'HGVSp',
  labelPos =c(406,428,515,351,356,361),
  legendTxtSize=1.2,
  showDomainLabel=TRUE
)

###################################################################################################
######  Gemcitabine nab paclitaxel response related major gene : TP53 lollipopPlot  ###############
###################################################################################################
AG_MAFs<- merge_mafs(c(rAG_mafs,mrAG_mafs, nrAG_mafs))
lollipopPlot(
  maf = AG_MAFs,
  gene = 'TP53',
  AACol= 'HGVSp',
  #labelPos =c(406,428,515,351,356,361),
  legendTxtSize=1.2,
  showDomainLabel=TRUE
)

###################################################################################################
#################  SMAD4: Association between LGD mutation domains and PFS  #######################
###################################################################################################
meta <- meta_YUHS
PDAC_lib_df <- PDAC_lib_df_YUHS
PDAC_MAFs <- PDAC_MAFs_YUHS
PDAC_VAF_mat <- PDAC_VAF_mat_YUHS

is_outlier <- function(x) {
  return(x < quantile(x, 0.4) - 0.5 * IQR(x) | x > quantile(x, 0.6) + 0.5 * IQR(x))
}

check_meta <- meta
check_meta$group <- 'wild-type' ; check_meta$group[check_meta$SMAD4==1] <-'non-LGD'
check_meta$group[check_meta$SAMPLE]
LGD_nonDomain1=sapply(strsplit(unique(as.vector(PDAC_MAFs@data[PDAC_MAFs@data$Hugo_Symbol %in% 'SMAD4' & (PDAC_MAFs@data$Variant_Classification %in% c('Frame_Shift_Del','Frame_Shift_Ins', 'Splice_Site','Nonsense_Mutation')),'Source_MAF' ])$Source_MAF ), "\\."), function(x) x[[1]])
LGD_nondomain=c()
for( i in seq(length(LGD_nonDomain1))){
  sam=LGD_nonDomain1[i]
  if(grepl("^[C]", sam)){
    LGD_nondomain <- c(LGD_nondomain,strsplit(unique(sam),'_')[[1]][1])
  }else if(grepl("^[P_]", sam)){
    LGD_nondomain <- c(LGD_nondomain,strsplit(unique(sam),'-')[[1]][1])
  }else{
    LGD_nondomain <- c(LGD_nondomain,strsplit(sam,'D')[[1]][1])
  }
}; LGD_nondomain

LGD_tot=sapply(strsplit(unique(as.vector(PDAC_MAFs@data[PDAC_MAFs@data$Hugo_Symbol %in% 'SMAD4' & PDAC_MAFs@data$Protein_position<=142 & PDAC_MAFs@data$Protein_position>=1 & (PDAC_MAFs@data$Variant_Classification %in% c('Frame_Shift_Del','Frame_Shift_Ins', 'Splice_Site','Nonsense_Mutation')),'Source_MAF' ])$Source_MAF ), "\\."), function(x) x[[1]])
LGD=c()
for( i in seq(length(LGD_tot))){
  sam=LGD_tot[i]
  if(grepl("^[C]", sam)){
    LGD <- c(LGD,strsplit(unique(sam),'_')[[1]][1])
  }else if(grepl("^[P_]", sam)){
    LGD <- c(LGD,strsplit(unique(sam),'-')[[1]][1])
  }else{
    LGD <- c(LGD,strsplit(sam,'D')[[1]][1])
  }
}; LGD
check_meta$group[paste0(check_meta$SAMPLE) %in% LGD_nondomain] <-'LGD+non-DBD'
check_meta$group[paste0(check_meta$SAMPLE) %in% LGD] <-'LGD+DBD'
check_meta <- check_meta[!is.na(check_meta$FFX_PFS),]
check_meta_1 <- check_meta[check_meta$group=='LGD+DBD',] ; check_meta_1 <- check_meta_1[!(is_outlier(check_meta_1$FFX_PFS)),]
check_meta_2 <- check_meta[check_meta$group=='LGD+non-DBD',] ; check_meta_2 <- check_meta_2[!(is_outlier(check_meta_2$FFX_PFS)),]
check_meta_3 <- check_meta[check_meta$group=='non-LGD',] ; check_meta_3 <- check_meta_3[!(is_outlier(check_meta_3$FFX_PFS)),]
check_meta_4 <- check_meta[check_meta$group=='wild-type',] ; check_meta_4 <- check_meta_4[!(is_outlier(check_meta_4$FFX_PFS)),]
check_meta <- rbind(check_meta_1,check_meta_2,check_meta_3) ; check_meta <- check_meta[!is.na(check_meta$FFX_PFS),]
check_meta <- check_meta[check_meta$group %in% c('LGD+DBD','LGD+non-DBD','non-LGD'),]
check_meta$group <- factor(check_meta$group)
fit<-survfit(Surv(FFX_PFS)~group, data=check_meta)
YUHS_res3=survdiff(Surv(FFX_PFS)~group, data=check_meta)$p
ggplot(check_meta, aes(x = factor(group), y = FFX_PFS, fill = factor(group))) +
  geom_boxplot() +
  labs(title = "FFX_PFS by Group", x = "Group", y = "FFX_PFS (months)") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("LGD+DBD" = "#F45050", "LGD+non-DBD" = "#B1AFFF", "non-LGD" = "#7C73C0", 'wild-type'='gray'))+
  geom_dotplot(method   = "dotdensity",
               alpha    = 1,
               binwidth = (max(check_meta$FFX_PFS) - min(check_meta$FFX_PFS))/70,
               colour   = "black",
               binaxis  = 'y',
               stackdir = 'center') +
  labs(x     = "Groups",
       y     = "PFS (days)") +
  stat_compare_means(aes(group = group), 
                     size = 5) +
  theme_bw() +
  theme(axis.title       = element_text(size = 25,
                                        face = 'plain'),
        title            = element_text(size = 25,
                                        face = 'italic'),
        plot.title       = element_blank(),
        axis.ticks = element_line(size = 1,color = 'black'),
        axis.ticks.length = unit(.2, "cm"),
        plot.background  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "black",
                                        fill   = NA,
                                        size   = 1,),
        legend.position = "bottom",
        legend.title = element_blank()) +
  theme(axis.line = element_line(color = 'black'),
        axis.text.x = element_text(size = 25,
                                   color = "black"),
        axis.text.y = element_text(size = 25,
                                   color = "black"))

###################################################################################################
#################  TP53: Association between LGD mutation domains and PFS  #######################
###################################################################################################
meta <- meta_YUHS
PDAC_lib_df <- PDAC_lib_df_YUHS
PDAC_MAFs <- PDAC_MAFs_YUHS
PDAC_VAF_mat <- PDAC_VAF_mat_YUHS
is_outlier <- function(x) {
  return(x < quantile(x, 0.4) - 0.5 * IQR(x) | x > quantile(x, 0.6) + 0.5 * IQR(x))
}
check_meta <- meta
check_meta$group <- 'wild-type' ; check_meta$group[check_meta$TP53==1] <-'non-LGD'
check_meta$group[check_meta$SAMPLE]
LGD_nonDomain1=sapply(strsplit(unique(as.vector(PDAC_MAFs@data[PDAC_MAFs@data$Hugo_Symbol %in% 'TP53' & (PDAC_MAFs@data$Variant_Classification %in% c('Frame_Shift_Del','Frame_Shift_Ins', 'Splice_Site','Nonsense_Mutation')),'Source_MAF' ])$Source_MAF ), "\\."), function(x) x[[1]])
LGD_nondomain=c()
for( i in seq(length(LGD_nonDomain1))){
  sam=LGD_nonDomain1[i]
  if(grepl("^[C]", sam)){
    LGD_nondomain <- c(LGD_nondomain,strsplit(unique(sam),'_')[[1]][1])
  }else if(grepl("^[P_]", sam)){
    LGD_nondomain <- c(LGD_nondomain,strsplit(unique(sam),'-')[[1]][1])
  }else{
    LGD_nondomain <- c(LGD_nondomain,strsplit(sam,'D')[[1]][1])
  }
}; LGD_nondomain

LGD_tot=sapply(strsplit(unique(as.vector(PDAC_MAFs@data[PDAC_MAFs@data$Hugo_Symbol %in% 'TP53' &PDAC_MAFs@data$Protein_position>=100 & PDAC_MAFs@data$Protein_position<=292 & (PDAC_MAFs@data$Variant_Classification %in% c('Frame_Shift_Del','Frame_Shift_Ins', 'Splice_Site','Nonsense_Mutation')),'Source_MAF' ])$Source_MAF ), "\\."), function(x) x[[1]])
LGD=c()
for( i in seq(length(LGD_tot))){
  sam=LGD_tot[i]
  if(grepl("^[C]", sam)){
    LGD <- c(LGD,strsplit(unique(sam),'_')[[1]][1])
  }else if(grepl("^[P_]", sam)){
    LGD <- c(LGD,strsplit(unique(sam),'-')[[1]][1])
  }else{
    LGD <- c(LGD,strsplit(sam,'D')[[1]][1])
  }
}; LGD
check_meta$group[paste0(check_meta$SAMPLE) %in% LGD_nondomain] <-'LGD+non-DBD'
check_meta$group[paste0(check_meta$SAMPLE) %in% LGD] <-'LGD+DBD'
check_meta <- check_meta[!is.na(check_meta$AG_PFS),]
check_meta_1 <- check_meta[check_meta$group=='LGD+DBD',] ; check_meta_1 <- check_meta_1[!(is_outlier(check_meta_1$AG_PFS)),]
check_meta_2 <- check_meta[check_meta$group=='LGD+non-DBD',] ; check_meta_2 <- check_meta_2[!(is_outlier(check_meta_2$AG_PFS)),]
check_meta_3 <- check_meta[check_meta$group=='non-LGD',] ; check_meta_3 <- check_meta_3[!(is_outlier(check_meta_3$AG_PFS)),]
check_meta_4 <- check_meta[check_meta$group=='wild-type',] ; check_meta_4 <- check_meta_4[!(is_outlier(check_meta_4$AG_PFS)),]
check_meta <- rbind(check_meta_1,check_meta_2,check_meta_3) ; check_meta <- check_meta[!is.na(check_meta$AG_PFS),]
check_meta <- check_meta[check_meta$group %in% c('LGD+DBD','LGD+non-DBD','non-LGD','wild-type'),]
check_meta$group <- factor(check_meta$group)
fit<-survfit(Surv(AG_PFS)~group, data=check_meta)
YUHS_res3=survdiff(Surv(AG_PFS)~group, data=check_meta)$p
ggplot(check_meta, aes(x = factor(group), y = AG_PFS, fill = factor(group))) +
  geom_boxplot() +
  labs(title = "FFX_PFS by Group", x = "Group", y = "FFX_PFS (months)") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("LGD+DBD" = "#F45050", "LGD+non-DBD" = "#B1AFFF", "non-LGD" = "#7C73C0", 'wild-type'='gray'))+
  geom_dotplot(method   = "dotdensity",
               alpha    = 1,
               binwidth = (max(check_meta$AG_PFS) - min(check_meta$AG_PFS))/70,
               colour   = "black",
               binaxis  = 'y',
               stackdir = 'center') +
  labs(x     = "Groups",
       y     = "PFS (days)") +
  stat_compare_means(aes(group = group), 
                     size = 5) +
  theme_bw() +
  theme(axis.title       = element_text(size = 25,
                                        face = 'plain'),
        title            = element_text(size = 25,
                                        face = 'italic'),
        plot.title       = element_blank(),
        axis.ticks = element_line(size = 1,color = 'black'),
        axis.ticks.length = unit(.2, "cm"),
        plot.background  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "black",
                                        fill   = NA,
                                        size   = 1,),
        legend.position = "bottom",
        legend.title = element_blank()) +
  theme(axis.line = element_line(color = 'black'),
        axis.text.x = element_text(size = 25,
                                   color = "black"),
        axis.text.y = element_text(size = 25,
                                   color = "black"))

