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

set.seed(123)

setwd('/Users/kamran/Documents/Spring 2022/Research/Models/Week 4-7/All')

binary_table <- read.csv('/Users/kamran/Documents/Spring 2022/Research/Models/Week 3-17/binary_table.csv')

`%!in%` <- Negate(`%in%`)

binary_no_top <- binary_table[(binary_table$Mega_peaks == 0 & binary_table$MEP_peaks == 0 & binary_table$Ery_peaks ==0),]
binary_top <- binary_table[binary_table$snpID %!in% binary_no_top$snpID,]

#check
(nrow(binary_no_top) + nrow(binary_top) - nrow(binary_table)) == 0

get_auc <- function(table, string_name) {
  
  z <- table %>% select(-c('snpID','chr','pos', 'rsID','Mega_peaks','MEP_peaks','Ery_peaks'))
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
  
  #return AUC
  write.csv(summary(gbm_fit, n.trees = optimal_num_trees, plotit = FALSE), paste0("boost_perf",string_name,'.csv'), row.names = F)
  return(auc(roc(GWAS_test,gbm_probabilities)))
}

print(get_auc(binary_no_top,'none_drop'))
print(get_auc(binary_top,'some_drop'))