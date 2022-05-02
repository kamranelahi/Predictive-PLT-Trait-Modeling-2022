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

binary_table <- read.csv('../binary_table.csv', header = T)
blood_traits <- names(binary_table[,c(9:15,237:254)])

first_cols <- binary_table[,c(1:8)]

AUCs_boosting <- seq(.1,1,.1)

set.seed(123)

test_limit <- function(table, cutoff){
  
  #add interaction terms
  table_sub <- table[,names(table) %in% blood_traits]
  cor_binary <- cor(table_sub)
  n <- data.frame(table[,1])
  for (k in 1:ncol(cor_binary)){
    for (j in 1:ncol(cor_binary)){
      if (cor_binary[k,j] > cutoff && j!=k){
        n <- cbind(n, table[,colnames(cor_binary)[k]]*table[,colnames(cor_binary)[j]])
        colnames(n)[ncol(n)] <- paste0(colnames(cor_binary)[k],"*",colnames(cor_binary)[j])
      }
    }
  }
  
  #create table
  n <- n[,-1]
  z <- cbind(table,n) %>% select(-c('snpID','chr','pos', 'rsID'))
  shuffled_data= z[sample(1:nrow(z)), ]
  divider <- as.integer(0.8*nrow(shuffled_data))
  train_factors <- shuffled_data[1:divider,]
  test_factors <- shuffled_data[(divider+1):nrow(shuffled_data),]
  GWAS_test <- test_factors$GWAS
  test_factors$GWAS <- NULL
  
  #create model
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
  
  #return AUC
  write.csv(summary(gbm_fit, n.trees = optimal_num_trees, plotit = FALSE), paste0("./boost_perf", toString(cutoff), ".csv"), row.names = F)
  return(auc(roc(GWAS_test,gbm_probabilities)))
}

for (i in seq(.1,1,.1)){
  AUCs_boosting[i*10] <- test_limit(binary_table, i)
  print(i)
}

write.csv(AUCs_boosting,'AUCs.csv',row.names = F)
