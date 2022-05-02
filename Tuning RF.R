library(dplyr)
library(tidyverse)
library(pROC)
library(kernlab)

set.seed(100)

binary_table <- read.csv('/Users/kamran/Documents/Spring 2022/Research/Models/Week 3-17/binary_table.csv', header = TRUE)

z <- binary_table %>% select(-c('snpID','chr','pos', 'rsID'))
shuffled_data= z[sample(1:nrow(z)), ]
divider <- as.integer(0.8*nrow(shuffled_data))
train_factors <- shuffled_data[1:divider,]
test_factors <- shuffled_data[(divider+1):nrow(shuffled_data),]
GWAS_test <- test_factors$GWAS
test_factors$GWAS <- NULL

library(ranger)

hyperparameter_grid = expand.grid(
  mtry  = c(1,2,3,4,5),
  num.trees = seq(300,1500,200),
  min.node.size   = seq(1,19,3),
  replace         = c(TRUE, FALSE), 
  sample.fraction = c(0.20, 0.30, 0.5, 0.63, 0.8, 1.0), 
  AUCs           = NA
)

system.time({
  
  for(i in 1:nrow(hyperparameter_grid)) 
  {
    
    tmp = ranger(
      formula         = GWAS  ~ ., 
      data            = train_factors, 
      num.trees       = hyperparameter_grid$num.trees[i],     
      mtry            = hyperparameter_grid$mtry[i],
      min.node.size   = hyperparameter_grid$min.node.size[i],
      replace         = hyperparameter_grid$replace[i],
      sample.fraction = hyperparameter_grid$sample.fraction[i],
      verbose         = FALSE,
      probability = TRUE,
      seed            = 100,
    )
    
    y_hat = predict(tmp, data = test_factors)$predictions
    this_auc = auc(roc(GWAS_test,y_hat[,2]))
    hyperparameter_grid$AUCs[i] = this_auc
    print(paste(i,this_auc,sep=' '))
  }
  
})

best_combo_index = which.max(hyperparameter_grid$AUCs)
hyperparameter_grid[best_combo_index,]

system.time({
  final_model = ranger(
    formula         = GWAS ~ ., 
    data            = train_factors, 
    num.trees       = hyperparameter_grid[best_combo_index, "num.trees"],
    mtry            = hyperparameter_grid[best_combo_index, "mtry"],
    min.node.size   = hyperparameter_grid[best_combo_index, "min.node.size"],
    replace         = hyperparameter_grid[best_combo_index, "replace"],
    sample.fraction = hyperparameter_grid[best_combo_index, "sample.fraction"],
    verbose         = FALSE,
    importance = "permutation",
    seed            = 777 )
})

sorted_importance = sort(final_model$variable.importance, decreasing = TRUE)
imp = data.frame(sorted_importance)

write.csv(imp,'rf_imp.csv',row.names=TRUE)
