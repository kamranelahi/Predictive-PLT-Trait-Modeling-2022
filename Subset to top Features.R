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

binary_table <- read.csv('../binary_table.csv', header = T)
factors_to_sub <- read.csv('/boost_perf.csv', header = T)

for (i in 0:10){
  limit <- i*0.25
  factors_tmp <- factors_to_sub[factors_to_sub$rel.inf >= limit,]
  z <- binary_table[,names(binary_table) %in% factors_tmp$var]
  z <- cbind(binary_table[,c(1:8)], z)
  
  x <- data.matrix(z[,c(5:7,9:ncol(z))])
  y <- data.matrix(z[,8])
  n <- ncol(x)
  pfactor = matrix(nrow = 1, ncol = n)
  pfactor[1,] = 1
  pfactor[1,c(1:3)] = 0
  cvfit = cv.glmnet(x, y, penalty.factor = pfactor, family = "binomial", type.measure = "auc")
  
  pdf(paste0("subset",toString(limit), "LASSO.output.AUCplot.pdf"))
  plot(cvfit, main = toString(limit))
  dev.off()
  
  cat("print(cvfit)\n")
  print(cvfit)
  cat("coef(cvfit, s=lambda.1se)\n")
  coef(cvfit, s="lambda.1se")
  cat("coef(cvfit, s=lambda.min)\n")
  coef(cvfit, s="lambda.min")
}