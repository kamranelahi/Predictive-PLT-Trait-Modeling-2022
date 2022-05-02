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
`%!in%` <- Negate(`%in%`)

set.seed(123)

binary_table <- read.csv('./binary_table.csv', header = T)
blood_traits <- names(binary_table[,237:254])

first_cols <- binary_table[,c(1:8)]

AUCs <- data.frame(matrix(NA, nrow = 10, ncol = 11))
names(AUCs) <- seq(0,2.5,.25)
rownames(AUCs) <- seq(.1,1,.1)

#function for adding interaction terms
add_interactions <- function(cutoff){
  cor_binary <- cor(binary_table[,237:254])
  n <- data.frame(binary_table[,1])
  for (k in 1:ncol(cor_binary)){
    for (j in 1:ncol(cor_binary)){
      if (cor_binary[k,j] > cutoff && j!=k){
        n <- cbind(n, binary_table[,colnames(cor_binary)[k]]*binary_table[,colnames(cor_binary)[j]])
        colnames(n)[ncol(n)] <- paste0(colnames(cor_binary)[k],"*",colnames(cor_binary)[j])
      }
    }
  }
  n <- n[,-1]
  return(cbind(binary_table,n))
}

#function for getting boost grid
get_top_factors <- function(table) {
  no_descript <- table %>% select(-c('snpID','chr','pos','rsID'))
  gbm_fit = gbm(GWAS ~ .,
                  distribution = "bernoulli",
                  n.trees = 1000,
                  interaction.depth = 2,
                  shrinkage = 0.01,
                  cv.folds = 5,
                  data = no_descript)
  gbm_fit_optimal = gbm_fit
  optimal_num_trees = gbm.perf(gbm_fit, plot.it = FALSE)
  return(summary(gbm_fit_optimal, n.trees = optimal_num_trees, plotit = FALSE))
}

#run the actual models
for (l in seq(.1,1,.1)){
  table_w_interactions <- add_interactions(l)
  factors_scored <- get_top_factors(table_w_interactions)
  write.csv(factors_scored, paste0('./Boosting Factors/boosting_factors', l, '.csv'), row.names = F)
  
  for (i in 0:10){
    limit <- i*0.25
    top_factors <- factors_scored[factors_scored$rel.inf >= limit,]
    z <- table_w_interactions[,names(table_w_interactions) %in% top_factors$var]
    added <- first_cols[,names(first_cols) %!in% names(z)]
    z <- cbind(added, z)
    
    non_gen_factors <- cbind(z$snp_maf,z$dist_nearest_gene,z$friends_ld07)
    x <- cbind(non_gen_factors,data.matrix(z[,c(9:ncol(z))]))
    y <- data.matrix(z$GWAS)
    xcol <- ncol(x)
    pfactor = matrix(nrow = 1, ncol = xcol)
    pfactor[1,] = 1
    pfactor[1,c(1:3)] = 0
    cvfit = cv.glmnet(x, y, penalty.factor = pfactor, family = "binomial", type.measure = "auc")
    
    AUCs[10*l,i+1] <- max(cvfit$cvm)
    print(paste(i,l,sep=" "))
  }
}


