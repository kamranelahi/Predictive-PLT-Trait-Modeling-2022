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
get_top_factors <- function(table, i) {
  no_descript <- table %>% select(-c('snpID','chr','pos','rsID'))
  final_model = ranger(
    formula         = GWAS ~ ., 
    data            = no_descript, 
    num.trees       = 300,
    mtry            = 5,
    min.node.size   = 1,
    replace         = F,
    sample.fraction = 0.63,
    verbose         = FALSE,
    importance = "permutation",
    seed            = 777 )
  return(data.frame(sort(final_model$variable.importance, decreasing = TRUE)))
}

#run the actual models
for (l in seq(.1,1,.1)){
  table_w_interactions <- add_interactions(l)
  factors_scored <- get_top_factors(table_w_interactions)
  factors_scored[,2] <- factors_scored[,1]
  factors_scored[,1] <- row.names(factors_scored)
  write.csv(factors_scored, paste0('./RF_factors', l, '.csv'), row.names = F)
  
  for (i in 0:10){
    
    decile_size <- as.integer(nrow(factors_scored)/10)
    top_factors <- factors_scored[1:i*decile_size,]
    
    z <- table_w_interactions[,names(table_w_interactions) %in% top_factors[,1]]
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


