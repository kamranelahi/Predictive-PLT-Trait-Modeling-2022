library(dplyr)
library(tidyr)
library(glmnet)
library(ISLR)
library(gbm)
library(tidyverse)
library(randomForest)
library(adabag)
library(data.table)
library(xgboost)
library(pROC)

binary_table <- fread('/Users/kamran/Documents/Spring 2022/Research/Models/Week 3-17/binary_table.csv')
plt <- fread('/Users/kamran/Documents/Spring 2022/Research/SNP GWAS Data/RBC WBC/OneDrive_1_3-15-2022/PLT.quant.processed.bed', 
             select = c('name','pval'))
rbc <- fread('/Users/kamran/Documents/Spring 2022/Research/SNP GWAS Data/RBC WBC/OneDrive_1_3-15-2022/RBC.quant.processed.bed', 
             select = c('name','pval'))
wbc <- fread('/Users/kamran/Documents/Spring 2022/Research/SNP GWAS Data/RBC WBC/OneDrive_1_3-15-2022/WBC.quant.processed.bed', 
             select = c('name','pval'))
inputs <- fread('/Users/kamran/Documents/Spring 2022/Research/SNP GWAS Data/SNPsnap_Matched_data/input_snps_annotated.txt', 
)
matches <- fread('/Users/kamran/Documents/Spring 2022/Research/SNP GWAS Data/SNPsnap_Matched_data/matched_snps_annotated.txt',
)
matches$input_rsID <- 0

set.seed(123)


`%!in%` <- Negate(`%in%`)

for (i in 1:nrow(matches)){
  in_sub <- which(inputs$snpID == matches$input_snp[i], arr.ind = FALSE)
  matches$input_rsID[i] <- inputs$rsID[in_sub[1]]
}

plt <- plt[plt$name %in% inputs$rsID,]
rbc <- rbc[rbc$name %in% inputs$rsID,]
wbc <- wbc[wbc$name %in% inputs$rsID,]

#switch is trait that selects for trait being included
# 0 = PLT
# 1 = RBC
# 2 = WBC

remove_snps <- function(subscrpt, iter, trait){
  pval_cutoff <- 5/(10^subscrpt)
  
  r_sub <- rbc[rbc$pval < pval_cutoff,]
  w_sub <- wbc[wbc$pval < pval_cutoff,]
  p_sub <- plt[plt$pval < pval_cutoff,]
  
  exclusive_matches <- matches[matches$input_rsID %in% r_sub$name | (matches$input_rsID %in% w_sub$name | matches$input_rsID %in% p_sub$name),]
  
  toInclude <- append(append(r_sub$name, append(w_sub$name, p_sub$name)), exclusive_matches$rsID)
  
  binary_subbed <- binary_table[binary_table$rsID %in% toInclude,]
  
  print(paste(pval_cutoff, length(binary_table$GWAS) - length(binary_subbed$GWAS), sep = ' '))
  
  z <- binary_subbed %>% select(-c('snpID','chr','pos', 'rsID'))
  shuffled_data= z[sample(1:nrow(z)), ]
  divider <- as.integer(0.8*nrow(shuffled_data))
  train_factors <- shuffled_data[1:divider,]
  test_factors <- shuffled_data[(divider+1):nrow(shuffled_data),]
  GWAS_test <- test_factors$GWAS
  test_factors$GWAS <- NULL
  
  gbm_fit = gbm(GWAS ~ .,
                distribution = "bernoulli",
                n.trees = 1000,
                interaction.depth = 2,
                shrinkage = 0.01,
                cv.folds = 5,
                data = train_factors)
  
  optimal_num_trees = gbm.perf(gbm_fit, plot.it = FALSE)
  gbm_probabilities = predict(gbm_fit, n.trees = optimal_num_trees,
                              type = "response", newdata = test_factors)
  
  #get factors
  auc_here <- round(auc(roc(GWAS_test,gbm_probabilities)), digits = 3)
  write.csv(summary(gbm_fit, n.trees = optimal_num_trees, plotit = FALSE),paste0('./factors',subscrpt,auc_here,'.csv'), row.names = F)
  
  #return AUC
  return(auc_here)
  
}

getAUCs <- function(t){
  AUCs <- data.frame(matrix(NA, nrow = 10, ncol = 10))
  for (i in c(2:10)){
    for (j in 1:10){
      AUCs[j,i-1] <- remove_snps(i,j,t)
    }
    write.csv(data.frame(AUCs), paste0('./AUCs-',t,'.csv'), row.names = FALSE)
  }
  write.csv(data.frame(AUCs), paste0('./AUCs-',t,'.csv'), row.names = FALSE)
}

for (k in 0){
  getAUCs(k)
}