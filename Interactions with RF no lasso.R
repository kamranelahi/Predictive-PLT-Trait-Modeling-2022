library(dplyr)
library(tidyr)
library(glmnet)
library(ISLR)
library(gbm)
library(tidyverse)
library(adabag)
library(data.table)
library(xgboost)
library(pROC)
library(ranger)
`%!in%` <- Negate(`%in%`)

set.seed(123)

binary_table <- read.csv('/Users/kamran/Documents/Spring 2022/Research/Models/Week 3-17/binary_table.csv', header = T)
blood_traits <- names(binary_table[,237:254])

first_cols <- binary_table[,c(1:8)]

AUCs <- data.frame(matrix(NA, nrow = 10, ncol = 11))
names(AUCs) <- seq(0,2.5,.25)
rownames(AUCs) <- seq(.1,1,.1)

rf_aucs <- seq(1,10,1)

#function for adding interaction terms
add_interactions <- function(cutoff){
  cor_binary <- cor(binary_table[,237:254])
  n <- data.frame(binary_table[,1])
  for (k in 1:ncol(cor_binary)){
    for (j in 1:ncol(cor_binary)){
      if (cor_binary[k,j] > cutoff && j!=k){
        n <- cbind(n, binary_table[,colnames(cor_binary)[k]]*binary_table[,colnames(cor_binary)[j]])
        colnames(n)[ncol(n)] <- paste0(colnames(cor_binary)[k],"x",colnames(cor_binary)[j])
      }
    }
  }
  n <- n[,-1]
  return(cbind(binary_table,n))
}

#function for getting boost grid
get_top_factors <- function(table) {
  
  z <- table %>% select(-c('snpID','chr','pos', 'rsID'))
  shuffled_data= z[sample(1:nrow(z)), ]
  divider <- as.integer(0.8*nrow(shuffled_data))
  train_factors <- shuffled_data[1:divider,]
  test_factors <- shuffled_data[(divider+1):nrow(shuffled_data),]
  GWAS_test <- test_factors$GWAS
  test_factors$GWAS <- NULL
  
  final_model = ranger(
    formula         = GWAS ~ ., 
    data            = train_factors, 
    num.trees       = 300,
    mtry            = 5,
    min.node.size   = 1,
    replace         = F,
    sample.fraction = 0.63,
    verbose         = FALSE,
    seed            = 777 )
  
  y_hat = predict(final_model, data = test_factors)$predictions
  this_auc = auc(roc(GWAS_test,y_hat))
  
  return(this_auc)
}

#run the actual models
for (l in seq(.1,1,.1)){
  table_w_interactions <- add_interactions(l)
  rf <- get_top_factors(table_w_interactions)
  rf_aucs[l*10] <- rf
    print(paste(l))
}

