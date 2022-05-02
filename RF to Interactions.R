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
factors_to_sub <- read.csv('./boost_perf.csv', header = T)
blood_traits <- names(binary_table[,237:254])

first_cols <- binary_table[,c(1:8)]

AUCs <- data.frame(matrix(NA, nrow = 10, ncol = 11))
names(AUCs) <- seq(0,2.5,.25)
rownames(AUCs) <- seq(.1,1,.1)

test_limit <- function(table, cutoff){
  table_sub <- table[,names(table) %in% blood_traits]
  cor_binary <- cor(table_sub)
  n <- data.frame(table[,1])
  for (k in 1:ncol(cor_binary)){
    for (j in 1:ncol(cor_binary)){
      if (cor_binary[k,j] > cutoff && j!=k){
        n <- cbind(n, table[,colnames(cor_binary)[k]]*table[,colnames(cor_binary)[j]])
        colnames(n)[ncol(n)] <- paste0(colnames(cor_binary)[k],"x",colnames(cor_binary)[j])
      }
    }
  }
  n <- n[,-1]
  z <- cbind(table,n)
  non_gen_factors <- cbind(z$snp_maf,z$dist_nearest_gene,z$friends_ld07)
  x_ <- cbind(non_gen_factors,data.matrix(z[,c(9:ncol(z))]))
  y <- data.matrix(z$GWAS)
  n <- ncol(x_)
  pfactor = matrix(nrow = 1, ncol = n)
  pfactor[1,] = 1
  pfactor[1,c(1:3)] = 0
  cvfit = cv.glmnet(x_, y, penalty.factor = pfactor, family = "binomial", type.measure = "auc")
  return(max(cvfit$cvm))
}

for (i in 0:10){
  limit <- i*0.25
  factors_tmp <- factors_to_sub[factors_to_sub$sorted_importance >= limit,]
  z_ <- binary_table[,names(binary_table) %in% factors_tmp[,1]]
  added <- first_cols[,names(first_cols) %!in% names(z_)]
  z_ <- cbind(added, z_)

  for (l in seq(.1,1,.1)) {
    AUCs[10*l,i+1] <- test_limit(z_,limit)
    print(paste(i,l, sep = ' '))
  }
}
