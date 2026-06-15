pkg_control(c('ComplexHeatmap','stringr',"maftools","tidyverse","ggplot2",'psych','rcompanion','mosaicData','ggcorrplot','ggplot2','msigdbr','annotate','org.Hs.eg.db',
              'factoextra','survminer','readxl','TiMEx','glmnet','vcdExtra','glmnet','ggmosaic','survival','circlize','tm','GSgalgoR','lsa','arkhe','wordcloud'))

##########################################################################################################
#############################  SMAD4: Mutually exclusive mutation pattern  ###############################
##########################################################################################################
ddr_gene <- read.table("/home/eunu0605/Project/PDAC/Library/GeneList/DDR_genelist.txt"); ddr_gene <- ddr_gene$V1
### 1) co-occurence ###
mut_dat <- data.frame(GENE1=c(),GENE2=c(),PAIR=c(),PVAL=c(),ODDS=c())
pdac_mut_mat_bin_test <- pdac_mut_mat_bin[,colSums(pdac_mut_mat_bin)>=0]
for(i in seq(ncol(pdac_mut_mat_bin_test))){
  gene1 <- colnames(pdac_mut_mat_bin_test)[i]
  for(j in seq(ncol(pdac_mut_mat_bin_test))){
    gene2 <- colnames(pdac_mut_mat_bin_test)[j]
    tmp_mat <- pdac_mut_mat_bin_test[,c(i,j)]
    tmp_table <- matrix( c(sum(rowSums(tmp_mat)==2),sum(tmp_mat[rowSums(tmp_mat)==1 & tmp_mat[,1]==1,1]),sum(tmp_mat[rowSums(tmp_mat)==1 & tmp_mat[,2]==1,2]),sum(rowSums(tmp_mat)==0)) , nrow = 2)
    if(tmp_table[1,1]==0 | tmp_table[1,2]==0 | tmp_table[2,1]==0 | tmp_table[2,2]==0 ){
      tmp_table[1,1] <-  tmp_table[1,1]+1 ;  tmp_table[1,2]<-  tmp_table[1,2]+1 ; tmp_table[2,1]<-  tmp_table[2,1]+1  ; tmp_table[2,2]<-  tmp_table[2,2]+1
    }
    fisher <- fisher.test(tmp_table) ; tmp_fisher <- fisher$p.value
    tmp_odds_ratio <- tmp_table[1,1]*tmp_table[2,2]/tmp_table[1,2]/tmp_table[2,1]
    mut_dat <- rbind(mut_dat,data.frame(GENE1=gene1,GENE2=gene2,PAIR=paste0(gene1,'+',gene2),PVAL=tmp_fisher,ODDS=tmp_odds_ratio))
    if(j==ncol(pdac_mut_mat_bin_test)){
      mut_dat$PVAL[mut_dat$GENE1==gene1] <- p.adjust(mut_dat$PVAL[mut_dat$GENE1==gene1],method='bonferroni')
    }
    print(c(i,j))
  }
};mut_dat <- mut_dat[mut_dat$GENE1 < mut_dat$GENE2, ]; mut_dat <- mut_dat[mut_dat$ODDS>1 & mut_dat$PVAL<0.2, ] ; rownames(mut_dat) <- seq(nrow(mut_dat))

### 2) mutually exclusive ###
pdac_mut_mat_bin_test <- pdac_mut_mat_bin[,]
TiMExx=TiMEx(pdac_mut_mat_bin_test, 0.5, 0.05, 0.05)

load('/home/eunu0605/Project/PDAC/Library/Clinic_info/mut_dat.rda')
TiMEx_gene=TiMExx$genesSignif[[2]]$fdr
SMAD_mut_gene=c(TiMEx_gene[TiMEx_gene[,1]=='SMAD4',2],TiMEx_gene[TiMEx_gene[,2]=='SMAD4',1])[!(c(TiMEx_gene[TiMEx_gene[,1]=='SMAD4',2],TiMEx_gene[TiMEx_gene[,2]=='SMAD4',1]) %in% ddr_gene)];SMAD_mut_gene
SMAD_mut_ddr_gene=c(TiMEx_gene[TiMEx_gene[,1]=='SMAD4',2],TiMEx_gene[TiMEx_gene[,2]=='SMAD4',1])[(c(TiMEx_gene[TiMEx_gene[,1]=='SMAD4',2],TiMEx_gene[TiMEx_gene[,2]=='SMAD4',1]) %in% ddr_gene)];SMAD_mut_ddr_gene
colSums(meta[,c(SMAD_mut_gene,'SMAD4')]); colSums(meta[,SMAD_mut_ddr_gene]) 
colSums(meta[meta$FFX %in% c(1,2,3),c(SMAD_mut_gene,'SMAD4')]); colSums(meta[meta$FFX %in% c(1,2,3),SMAD_mut_ddr_gene])
colSums(meta[meta$FFX %in% c(1,2),c(SMAD_mut_gene,'SMAD4')]); colSums(meta[meta$FFX %in% c(1,2),SMAD_mut_ddr_gene])
pdac_mut_mat_bin_test <- pdac_mut_mat_bin
for ( i  in seq(nrow(TiMEx_gene))){
  tmp_mat <- pdac_mut_mat_bin_test[,c(TiMEx_gene[i,1],TiMEx_gene[i,2])]
  tmp_table <- matrix( c(sum(rowSums(tmp_mat)==2),sum(tmp_mat[rowSums(tmp_mat)==1 & tmp_mat[,1]==1,1]),sum(tmp_mat[rowSums(tmp_mat)==1 & tmp_mat[,2]==1,2]),sum(rowSums(tmp_mat)==0)) , nrow = 2)
  if(tmp_table[1,1]==0 | tmp_table[1,2]==0 | tmp_table[2,1]==0 | tmp_table[2,2]==0 ){
    tmp_table[1,1] <-  tmp_table[1,1]+1 ;  tmp_table[1,2]<-  tmp_table[1,2]+1 ; tmp_table[2,1]<-  tmp_table[2,1]+1  ; tmp_table[2,2]<-  tmp_table[2,2]+1
  }
  tmp_odds_ratio <- tmp_table[1,1]*tmp_table[2,2]/tmp_table[1,2]/tmp_table[2,1]
  mut_dat <- rbind(mut_dat,
                   data.frame(GENE1=TiMEx_gene[i,1],GENE2=TiMEx_gene[i,2],PAIR=paste0(TiMEx_gene[i,1],'+',TiMEx_gene[i,2]),PVAL=TiMExx$pvals[[2]]$fdr[i],ODDS=tmp_odds_ratio))
  
}; rownames(mut_dat) <- seq(nrow(mut_dat))
### 3) mut p_val<0.2 ###
TiMExx2=TiMEx(pdac_mut_mat_bin_test, 0.5, 0.2, 0.2)
aa <- TiMExx2 ; rm(TiMExx2)
bb=aa$genesSignif[[2]]$fdr   
for ( i  in seq(nrow(bb))){
  tmp_mat <- pdac_mut_mat_bin_test[,c(bb[i,1],bb[i,2])]
  tmp_table <- matrix( c(sum(rowSums(tmp_mat)==2),sum(tmp_mat[rowSums(tmp_mat)==1 & tmp_mat[,1]==1,1]),sum(tmp_mat[rowSums(tmp_mat)==1 & tmp_mat[,2]==1,2]),sum(rowSums(tmp_mat)==0)) , nrow = 2)
  
  if(sum((mut_dat$GENE1==bb[i,1] & mut_dat$GENE2==bb[i,2] & mut_dat$ODDS<1) |(mut_dat$GENE2==bb[i,1] & mut_dat$GENE1==bb[i,2]& mut_dat$ODDS<1))>0 ){
    print(c(bb[i,1],bb[i,2]))
  }else{
    print(c(i,nrow(mut_dat[(mut_dat$GENE1=='TGFBR2' &  mut_dat$GENE2=='SMAD4') & mut_dat$PVAL<0.05,])))
    if(tmp_table[1,1]==0 | tmp_table[1,2]==0 | tmp_table[2,1]==0 | tmp_table[2,2]==0 ){
      tmp_table[1,1] <-  tmp_table[1,1]+1 ;  tmp_table[1,2]<-  tmp_table[1,2]+1 ; tmp_table[2,1]<-  tmp_table[2,1]+1  ; tmp_table[2,2]<-  tmp_table[2,2]+1
    }
    tmp_odds_ratio <- tmp_table[1,1]*tmp_table[2,2]/tmp_table[1,2]/tmp_table[2,1]
    if(paste0(bb[i,1],'+',bb[i,2]) %in% mut_dat$PAIR ){
      idx=mut_dat$PVAL %in% paste0(bb[i,1],'+',bb[i,2])
      mut_dat$PVAL[idx] <- aa$pvals[[2]]$fdr[i]
    }
    else if(paste0(bb[i,2],'+',bb[i,1]) %in% mut_dat$PAIR ){
      idx=mut_dat$PVAL %in% paste0(bb[i,2],'+',bb[i,1])
      mut_dat$PVAL[idx] <- aa$pvals[[2]]$fdr[i]
    }
    else {mut_dat <- rbind(mut_dat,data.frame(GENE1=bb[i,1],GENE2=bb[i,2],PAIR=paste0(bb[i,1],'+',bb[i,2]),PVAL=aa$pvals[[2]]$fdr[i],ODDS=tmp_odds_ratio))}
    
  }
}
hazard_meta <- meta[,c('SAMPLE','FFX_PFS','OS')]
pdac_mut_mat_bin_test <- pdac_mut_mat_bin[,]
pdac_mut_mat_bin_test[,'SMAD4'] <- 0 ; pdac_mut_mat_bin_test[paste0(meta$SAMPLE[meta$SMAD4type==1],'D'),'SMAD4'] <- 1
for(i in seq(nrow(mut_dat))){
  gene1=mut_dat$GENE1[i]
  gene2=mut_dat$GENE2[i]
  tmp_mat <- pdac_mut_mat_bin_test[,c(gene1,gene2)]
  if(mut_dat$ODDS[i]>1){
    hazard_meta[,paste0(gene1,'+',gene2)] <- as.numeric(rowSums(tmp_mat) %in% c(2))[match(rownames(pdac_mut_mat_bin_test),paste0(hazard_meta$SAMPLE,'D'))]
  }else{
    hazard_meta[,paste0(gene1,'+',gene2)] <- as.numeric(rowSums(tmp_mat) %in% c(1,2))[match(rownames(pdac_mut_mat_bin_test),paste0(hazard_meta$SAMPLE,'D'))]
  }
  print(paste0(i,':',ncol(hazard_meta)))
}

mut_dat$HAZARD <- 0 ; mut_dat$Hazard_pval <- 0 ; mut_dat$LOWER <- 0 ; mut_dat$UPPER <- 0  
for(i in seq(ncol(hazard_meta)-3)){
  hazard_meta_tmp=hazard_meta[,c(1,2,3,i+3)] ; colnames(hazard_meta_tmp)[4] <- 'FACTOR'
  hazard_meta_tmp$FACTOR <- as.factor(hazard_meta_tmp$FACTOR)
  if(length(unique(hazard_meta_tmp$FACTOR))==1){
    mut_dat$HAZARD[i] <- 1
    mut_dat$Hazard_pval[i] <-1
    mut_dat$LOWER[i] <-1
    mut_dat$UPPER[i] <-1
    
  }else {
    res_tmp<-coxph(Surv(FFX_PFS)~FACTOR,data= hazard_meta_tmp) 
    mut_dat$HAZARD[i] <- exp(res_tmp$coefficients)
    mut_dat$Hazard_pval[i] <-summary(res_tmp)$waldtest['pvalue']
    mut_dat$LOWER[i] <-summary(res_tmp,conf.int=0.90)$conf.int[3]
    mut_dat$UPPER[i] <-summary(res_tmp,conf.int=0.90)$conf.int[4] 
  }
  print(i)
}
FMGs=c('SMAD4',SMAD_mut_ddr_gene, SMAD_mut_gene)
sample_order <- character(0)
check_gene <- character(0)
for (letter in FMGs[!(FMGs %in% 'SETBP1')]) {
  if(length(check_gene)==0){
    samples <- meta$SAMPLE[meta[, letter] == 1 ]
    sample_order <- c(sample_order, samples)
  }else if(length(check_gene)==1){
    samples <- meta$SAMPLE[meta[, letter] == 1 & meta[, check_gene] ==0 ]
    sample_order <- c(sample_order, samples)
  }else{
    samples <- meta$SAMPLE[meta[, letter] == 1 & rowSums(meta[, check_gene]) ==0 ]
    sample_order <- c(sample_order, samples)
  }
  check_gene=c(check_gene,letter)
}
sample_order <- paste0(sample_order, 'D')
transposed_data_ordered <- t(pdac_mut_mat_bin)[FMGs, match(sample_order, colnames(t(pdac_mut_mat_bin)))]
colors <- c("#EEEEEE", "#F45050")
heatmap_object <- Heatmap(
  transposed_data_ordered,
  col = colors,
  column_order = sample_order,
  column_title = "Gene with Mutually Exclusive Mutation Pattern of SMAD4",
  show_column_names = FALSE,
  row_order = FMGs,
  row_names_gp = gpar(fontsize = 12, side = "left"), 
  column_names_gp = gpar(fontface = "bold"),
  column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  row_title_gp = gpar(fontsize = 20), 
  show_heatmap_legend = FALSE 
) ; draw(heatmap_object)


###################################################################################################
##########################  SMAD4: OR and HR relationshipt scatterplot  ###########################
###################################################################################################
mut_dat_test <- mut_dat[mut_dat$PVAL<0.05 & mut_dat$ODDS<1 &( mut_dat$GENE1=='SMAD4' |  mut_dat$GENE2=='SMAD4'),]
mut_dat_test$log_odds <- -log(mut_dat_test$ODDS)
mut_dat_test$log_hazard <- log(mut_dat_test$HAZARD)

plot <- ggplot(mut_dat_test, aes(x = log_odds, y = log_hazard)) +
  geom_point(color = "#F45050",size = 5, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#F45050") +
  theme_minimal() +
  geom_hline(yintercept = 0, color = "black") +
  ggrepel::geom_text_repel(data = as.data.frame(mut_dat_test),
                           ggplot2::aes(label = mut_dat_test[,'PAIR']),
                           size = 3) +
  labs(x     = "-log(Odds Ratio)",
       y     = "log(Hazard ratio)",
       title = "") +
  theme_bw() +
  theme(axis.title       = element_text(size = 20,
                                        face = 'plain'),
        title            = element_text(size = 20,
                                        face = 'bold'),
        plot.title       = element_text(hjust = 0.5),
        plot.background  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "black",
                                        fill   = NA,
                                        size   = 1,),
        legend.position = "bottom",
        legend.title = element_blank()) +
  theme(axis.line = element_line(color = 'black'),
        axis.text.x = element_text(size = 10,
                                   color = "black"),
        axis.text.y = element_text(size = 20,
                                   color = "black")) ; plot

##########################################################################################################
#############################  TP53: Mutually exclusive mutation pattern  ###############################
##########################################################################################################
### 1) co-occurence ###
mut_dat <- data.frame(GENE1=c(),GENE2=c(),PAIR=c(),PVAL=c(),ODDS=c())
pdac_mut_mat_bin_test <- pdac_mut_mat_bin[,colSums(pdac_mut_mat_bin)>=0]
for(i in seq(ncol(pdac_mut_mat_bin_test))){
  gene1 <- colnames(pdac_mut_mat_bin_test)[i]
  for(j in seq(ncol(pdac_mut_mat_bin_test))){
    gene2 <- colnames(pdac_mut_mat_bin_test)[j]
    tmp_mat <- pdac_mut_mat_bin_test[,c(i,j)]
    tmp_table <- matrix( c(sum(rowSums(tmp_mat)==2),sum(tmp_mat[rowSums(tmp_mat)==1 & tmp_mat[,1]==1,1]),sum(tmp_mat[rowSums(tmp_mat)==1 & tmp_mat[,2]==1,2]),sum(rowSums(tmp_mat)==0)) , nrow = 2)
    if(tmp_table[1,1]==0 | tmp_table[1,2]==0 | tmp_table[2,1]==0 | tmp_table[2,2]==0 ){
      tmp_table[1,1] <-  tmp_table[1,1]+1  ;  tmp_table[1,2]<-  tmp_table[1,2]+1 ; tmp_table[2,1]<-  tmp_table[2,1]+1  ; tmp_table[2,2]<-  tmp_table[2,2]+1
    }
    fisher <- fisher.test(tmp_table) ; tmp_fisher <- fisher$p.value
    tmp_odds_ratio <- tmp_table[1,1]*tmp_table[2,2]/tmp_table[1,2]/tmp_table[2,1]
    mut_dat <- rbind(mut_dat,data.frame(GENE1=gene1,GENE2=gene2,PAIR=paste0(gene1,'+',gene2),PVAL=tmp_fisher,ODDS=tmp_odds_ratio))
    if(j==ncol(pdac_mut_mat_bin_test)){
      mut_dat$PVAL[mut_dat$GENE1==gene1] <- p.adjust(mut_dat$PVAL[mut_dat$GENE1==gene1],method='bonferroni')
    }
    print(c(i,j))
  }
};mut_dat <- mut_dat[mut_dat$GENE1 < mut_dat$GENE2, ]; mut_dat <- mut_dat[mut_dat$ODDS>1 & mut_dat$PVAL<0.2, ] ; rownames(mut_dat) <- seq(nrow(mut_dat))

### 2) mut p_val<0.05 ###
pdac_mut_mat_bin_test <- pdac_mut_mat_bin[,colSums(pdac_mut_mat_bin)>=5]
TiMExx=TiMEx(pdac_mut_mat_bin_test, 0.5, 0.05, 0.05)

load('/home/eunu0605/Project/PDAC/Library/Clinic_info/mut_dat.rda')
TiMEx_gene=TiMExx$genesSignif[[2]]$fdr
TP53_mut_gene=c(TiMEx_gene[TiMEx_gene[,1]=='TP53',2],TiMEx_gene[TiMEx_gene[,2]=='TP53',1]);TP53_mut_gene
TP53_mut_gene=TP53_mut_gene[colSums(meta[,TP53_mut_gene])>=13]
colSums(meta[,c(TP53_mut_gene,'TP53')])
colSums(meta[meta$AG %in% c(1,2,3),c(TP53_mut_gene,'TP53')])
pdac_mut_mat_bin_test <- pdac_mut_mat_bin
for ( i  in seq(nrow(TiMEx_gene))){
  tmp_mat <- pdac_mut_mat_bin_test[,c(TiMEx_gene[i,1],TiMEx_gene[i,2])]
  tmp_table <- matrix( c(sum(rowSums(tmp_mat)==2),sum(tmp_mat[rowSums(tmp_mat)==1 & tmp_mat[,1]==1,1]),sum(tmp_mat[rowSums(tmp_mat)==1 & tmp_mat[,2]==1,2]),sum(rowSums(tmp_mat)==0)) , nrow = 2)
  if(tmp_table[1,1]==0 | tmp_table[1,2]==0 | tmp_table[2,1]==0 | tmp_table[2,2]==0 ){
    tmp_table[1,1] <-  tmp_table[1,1]+1 ;  tmp_table[1,2]<-  tmp_table[1,2]+1 ; tmp_table[2,1]<-  tmp_table[2,1]+1  ; tmp_table[2,2]<-  tmp_table[2,2]+1
  }
  tmp_odds_ratio <- tmp_table[1,1]*tmp_table[2,2]/tmp_table[1,2]/tmp_table[2,1]
  mut_dat <- rbind(mut_dat,data.frame(GENE1=TiMEx_gene[i,1],GENE2=TiMEx_gene[i,2],PAIR=paste0(TiMEx_gene[i,1],'+',TiMEx_gene[i,2]),PVAL=TiMExx$pvals[[2]]$fdr[i],ODDS=tmp_odds_ratio))
}; rownames(mut_dat) <- seq(nrow(mut_dat))

### 3) mut p_val<0.2 ###
TiMExx2=TiMEx(pdac_mut_mat_bin_test, 0.5, 0.2, 0.2)
aa <- TiMExx2 ; rm(TiMExx2)
bb=aa$genesSignif[[2]]$fdr  
for ( i  in seq(nrow(bb))){
  tmp_mat <- pdac_mut_mat_bin_test[,c(bb[i,1],bb[i,2])]
  tmp_table <- matrix( c(sum(rowSums(tmp_mat)==2),sum(tmp_mat[rowSums(tmp_mat)==1 & tmp_mat[,1]==1,1]),sum(tmp_mat[rowSums(tmp_mat)==1 & tmp_mat[,2]==1,2]),sum(rowSums(tmp_mat)==0)) , nrow = 2)
  
  if(sum((mut_dat$GENE1==bb[i,1] & mut_dat$GENE2==bb[i,2] & mut_dat$ODDS<1) |(mut_dat$GENE2==bb[i,1] & mut_dat$GENE1==bb[i,2]& mut_dat$ODDS<1))>0 ){
    print(c(bb[i,1],bb[i,2]))
  }else{
    print(c(i,nrow(mut_dat[(mut_dat$GENE1=='TGFBR2' &  mut_dat$GENE2=='SMAD4') & mut_dat$PVAL<0.05,])))
    if(tmp_table[1,1]==0 | tmp_table[1,2]==0 | tmp_table[2,1]==0 | tmp_table[2,2]==0 ){
      tmp_table[1,1] <-  tmp_table[1,1]+1 ;  tmp_table[1,2]<-  tmp_table[1,2]+1 ; tmp_table[2,1]<-  tmp_table[2,1]+1  ; tmp_table[2,2]<-  tmp_table[2,2]+1
    }
    tmp_odds_ratio <- tmp_table[1,1]*tmp_table[2,2]/tmp_table[1,2]/tmp_table[2,1]
    if(paste0(bb[i,1],'+',bb[i,2]) %in% mut_dat$PAIR ){
      idx=mut_dat$PVAL %in% paste0(bb[i,1],'+',bb[i,2])
      mut_dat$PVAL[idx] <- aa$pvals[[2]]$fdr[i]
    }
    else if(paste0(bb[i,2],'+',bb[i,1]) %in% mut_dat$PAIR ){
      idx=mut_dat$PVAL %in% paste0(bb[i,2],'+',bb[i,1])
      mut_dat$PVAL[idx] <- aa$pvals[[2]]$fdr[i]
    }
    else {mut_dat <- rbind(mut_dat,data.frame(GENE1=bb[i,1],GENE2=bb[i,2],PAIR=paste0(bb[i,1],'+',bb[i,2]),PVAL=aa$pvals[[2]]$fdr[i],ODDS=tmp_odds_ratio))}
    
  }
}
hazard_meta <- meta[meta$AG %in% c(1,3),c('SAMPLE','AG_PFS','OS')]
pdac_mut_mat_bin_test <- pdac_mut_mat_bin[paste0(c(rAG,nrAG),'D'),]

for(i in seq(nrow(mut_dat))){
  gene1=mut_dat$GENE1[i]
  gene2=mut_dat$GENE2[i]
  tmp_mat <- pdac_mut_mat_bin_test[,c(gene1,gene2)]
  if(mut_dat$ODDS[i]>1){
    hazard_meta[,paste0(gene1,'+',gene2)] <- as.numeric(rowSums(tmp_mat) %in% c(2))[match(rownames(pdac_mut_mat_bin_test),paste0(hazard_meta$SAMPLE,'D'))]
  }else{
    hazard_meta[,paste0(gene1,'+',gene2)] <- as.numeric(rowSums(tmp_mat) %in% c(1,2))[match(rownames(pdac_mut_mat_bin_test),paste0(hazard_meta$SAMPLE,'D'))]
  }
  print(paste0(i,':',ncol(hazard_meta)))
}
mut_dat$HAZARD <- 0 ; mut_dat$Hazard_pval <- 0 ; mut_dat$LOWER <- 0 ; mut_dat$UPPER <- 0  
for(i in seq(ncol(hazard_meta)-3)){
  hazard_meta_tmp=hazard_meta[,c(1,2,3,i+3)] ; colnames(hazard_meta_tmp)[4] <- 'FACTOR'
  hazard_meta_tmp$FACTOR <- as.factor(hazard_meta_tmp$FACTOR)
  if(length(unique(hazard_meta_tmp$FACTOR))==1){
    mut_dat$HAZARD[i] <- 1
    mut_dat$Hazard_pval[i] <-1
    mut_dat$LOWER[i] <-1
    mut_dat$UPPER[i] <-1
    
  }else {
    res_tmp<-coxph(Surv(AG_PFS)~FACTOR,data=hazard_meta_tmp) 
    mut_dat$HAZARD[i] <- exp(res_tmp$coefficients)
    mut_dat$Hazard_pval[i] <-summary(res_tmp)$waldtest['pvalue']
    mut_dat$LOWER[i] <-summary(res_tmp)$conf.int[3]
    mut_dat$UPPER[i] <-summary(res_tmp)$conf.int[4] 
  }
  print(i)
}
FMGs=c('TP53',TP53_mut_gene,'KRAS')
sample_order <- character(0)
check_gene <- character(0)
for (letter in FMGs[!(FMGs %in% 'SETBP1')]) {
  if(length(check_gene)==0){
    samples <- meta$SAMPLE[meta[, letter] == 1 ]
    sample_order <- c(sample_order, samples)
  }else if(length(check_gene)==1){
    samples <- meta$SAMPLE[meta[, letter] == 1 & meta[, check_gene] ==0 ]
    sample_order <- c(sample_order, samples)
  }else{
    samples <- meta$SAMPLE[meta[, letter] == 1 & rowSums(meta[, check_gene]) ==0 ]
    sample_order <- c(sample_order, samples)
  }
  check_gene=c(check_gene,letter)
}
sample_order <- paste0(sample_order, 'D')
transposed_data_ordered <- t(pdac_mut_mat_bin)[FMGs, match(sample_order, colnames(t(pdac_mut_mat_bin)))]
colors <- c("#EEEEEE", "#A8D1D1")
heatmap_object <- Heatmap(
  transposed_data_ordered,
  col = colors,
  column_order = sample_order,
  column_title = "Gene with mutually exclusive and co-occurence mutation pattern of TP53",
  show_column_names = FALSE,
  row_order = FMGs,
  row_names_gp = gpar(fontsize = 12, side = "left"), #
  column_names_gp = gpar(fontface = "bold"),
  column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  row_title_gp = gpar(fontsize = 20), 
  show_heatmap_legend = FALSE 
) ; draw(heatmap_object)


###################################################################################################
##########################  TP53: OR and HR relationshipt scatterplot  ###########################
###################################################################################################
mut_dat_test <- mut_dat[mut_dat$PVAL<0.05 & mut_dat$ODDS<1 &( (mut_dat$GENE1=='TP53' &mut_dat$GENE2 %in% TP53_mut_gene ) |  (mut_dat$GENE2=='TP53' &mut_dat$GENE1 %in% TP53_mut_gene )),]
mut_dat_test$log_odds <- -log(mut_dat_test$ODDS)
mut_dat_test$log_hazard <- log(mut_dat_test$HAZARD)

plot <- ggplot(mut_dat_test, aes(x = log_odds, y = log_hazard)) +
  geom_point(color = "#4AA96C",size = 5, alpha = 0.5) +
  theme_minimal() +
  geom_hline(yintercept = 0, color = "black") +
  ggrepel::geom_text_repel(data = as.data.frame(mut_dat_test),
                           ggplot2::aes(label = mut_dat_test[,'PAIR']),
                           size = 3) +
  labs(x     = "-log(Odds Ratio)",
       y     = "log(Hazard ratio)",
       title = "") +
  theme_bw() +
  theme(axis.title       = element_text(size = 20,
                                        face = 'plain'),
        title            = element_text(size = 20,
                                        face = 'bold'),
        plot.title       = element_text(hjust = 0.5),
        plot.background  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "black",
                                        fill   = NA,
                                        size   = 1,),
        legend.position = "bottom",
        legend.title = element_blank()) +
  theme(axis.line = element_line(color = 'black'),
        axis.text.x = element_text(size = 10,
                                   color = "black"),
        axis.text.y = element_text(size = 20,
                                   color = "black")) ; plot
