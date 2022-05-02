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

set.seed(100)
`%!in%` <- Negate(`%in%`)

setwd('/Users/kamran/Documents/Spring 2022/Research/Models/Week 4-7/')

binary_table <- read.csv('/Users/kamran/Documents/Spring 2022/Research/Models/Week 3-17/binary_table.csv')

FACTORS_TO_SPLIT = 25

get_boosting_model <- function(table) {
  
  z <- table %>% select(-c('snpID','chr','pos', 'rsID'))
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
  list_factors <- summary(gbm_fit, n.trees = optimal_num_trees, plotit = FALSE)
  list_factors <- rbind(list_factors, c('AUC',auc(roc(GWAS_test,gbm_probabilities))))
  return(list_factors)
  
}

original_list_factor <- get_boosting_model(binary_table)
yes <- c(original_list_factor[nrow(original_list_factor),2])
no <- c(original_list_factor[nrow(original_list_factor),2])
write.csv(original_list_factor,'all_pref.csv',row.names = F)

AUCs <- data.frame(yes,no)
rownames(AUCs) <- c('all')

top_factors <- (original_list_factor[original_list_factor[,1] %!in% c('snp_maf', 'dist_nearest_gene', 'friends_ld07', 'AUC'),])[1:FACTORS_TO_SPLIT,]

for (i in 1:nrow(top_factors)) {
  string_name <- top_factors[i,1]
  binary_no <- binary_table[binary_table[,string_name]==0,]
  binary_yes <- binary_table[binary_table[,string_name]==1,]
  list_no <- get_boosting_model(binary_no)
  list_yes <- get_boosting_model(binary_yes)[order(var),]
  write.csv(list_no,paste0(string_name,'_no.csv'),row.names = F)
  write.csv(list_yes,paste0(string_name,'_yes.csv'),row.names = F)
  difference <- data.frame(list_no$var, abs(list_no$rel.inf - list_yes$rel.inf))
  names(difference) <- c('var','abs_dif')
  write.csv(difference, paste0(string_name,'_difference.csv',row.names=F))
  new <- c(list_yes[nrow(list_yes),2],list_no[nrow(list_no),2])
  AUCs <- rbind(AUCs, new)
  rownames(AUCs)[i+1] <- string_name
}

write.csv(AUCs, 'AUCs.csv', row.names = TRUE)