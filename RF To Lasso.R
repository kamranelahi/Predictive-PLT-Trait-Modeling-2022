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

#binary_table <- read.csv('../binary_table.csv', header = T)
factors_to_sub <- read.csv('../../RF To Lasso/rf_imp.csv', header = T)


for (i in 1:10){
  
  decile_size <- as.integer(nrow(factors_to_sub)/10)
  
  factors_tmp <- factors_to_sub[1:i*decile_size,]
  z <- binary_table[,names(binary_table) %in% factors_tmp[,1]]
  z <- cbind(binary_table[,c(1:8)], z)
  
  x <- data.matrix(z[,c(5:7,9:ncol(z))])
  y <- data.matrix(z[,8])
  n <- ncol(x)
  pfactor = matrix(nrow = 1, ncol = n)
  pfactor[1,] = 1
  pfactor[1,c(1:3)] = 0
  cvfit = cv.glmnet(x, y, penalty.factor = pfactor, family = "binomial", type.measure = "auc")
  
  pdf(paste0("subset",toString(i), "LASSO.output.AUCplot.pdf"))
  plot(cvfit, main = toString(i))
  dev.off()
  
  cat("print(cvfit)\n")
  print(cvfit)
  cat("coef(cvfit, s=lambda.1se)\n")
  coef(cvfit, s="lambda.1se")
  cat("coef(cvfit, s=lambda.min)\n")
  coef(cvfit, s="lambda.min")
}