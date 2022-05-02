library(dplyr)
library(tidyr)
library(glmnet)
library(ISLR)
library(gbm)
library(tidyverse)
library(randomForest)
library(adabag)
library(data.table)s
library(xgboost)
library(pROC)
binary_table <- read.csv('/Users/kamran/Documents/Spring 2022/Research/Models/Week 3-17/binary_table.csv', header = T, stringsAsFactors = F)
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

AUCs <- data.frame(matrix(NA, nrow = 10, ncol = 10))

rbc <- rbc[rbc$name %in% inputs$rsID,]
wbc <- wbc[wbc$name %in% inputs$rsID,]

remove_wbc_rbc <- function(subscrpt, iter){
  pval_cutoff <- 5/(10^subscrpt)

  rbc_sub <- rbc[rbc$pval < pval_cutoff,]
  wbc_sub <- wbc[wbc$pval < pval_cutoff,]

  excluded_matches <- matches[matches$input_rsID %in% rbc_sub$name | matches$input_rsID %in% wbc_sub$name,]

  toRemove <- append(append(rbc_sub$name, wbc_sub$name), excluded_matches$rsID)

  binary_subbed <- binary_table[binary_table$rsID %!in% toRemove,]

  print(paste(pval_cutoff, length(binary_table[,1]) - length(binary_subbed[,1]), sep = ' '))
  
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
  if (iter == 1){
    write.csv(summary(gbm_fit, n.trees = optimal_num_trees, plotit = FALSE),paste0('./factors',subscrpt,'.csv'), row.names = F)
  } else{
    tmp <- read.csv(paste0('./factors',subscrpt,'.csv'))
    tmp <- rbind(tmp, summary(gbm_fit, n.trees = optimal_num_trees, plotit = FALSE))
    write.csv(tmp, paste0('./factors',subscrpt,'.csv'), row.names = F)
    print(nrow(tmp))
  }
  if (iter == 10){2
    write.csv(aggregate(tmp$rel.inf, list(tmp$var), FUN=mean), paste0('./factors',subscrpt,'.csv'), row.names = F)
  }
  
  #return AUC
  return(auc(roc(GWAS_test,gbm_probabilities)))
}

for (i in c(2:10,1000000)){
  for (j in 1:10){
    AUCs[j,i-1] <- remove_wbc_rbc(i, j)
  }
}

write.csv(data.frame(AUCs), './AUCs_no_wbc_rbc.csv', row.names = FALSE)
