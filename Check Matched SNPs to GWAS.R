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
matched_snps <- "matched_snps_annotated.txt"
gwas_table <- "PLT.Full.GWAS.bed"
p_value <- 0.05

filter_matched_snps <- function(matched_snps, gwas_table, p_value) {
  matched <- read.table(matched_snps, sep = '\t', header = T) %>% select('snpID','rsID')
  GWAS <- fread(gwas_table, select = c('name','pval'))
  GWAS_sub <- GWAS[GWAS$pval <= p_value,]
  matched$include <- ifelse(matched$rsID %in% GWAS_sub$name, 1,0)
  matched_exclude <- matched[matched$include == 1,]
  write.csv(matched_exclude, "excluded_matches.csv", row.names = F)
}

pvs <- c()
lgs <- c()
for (i in 1:10) {
  p_value <- p_value/10
  GWAS_sub <- GWAS_sub[GWAS_sub$pval <= p_value,]
  matched_exclude <- matched[matched$rsID %in% GWAS_sub$name,]
  pvs <- append(pvs, p_value)
  lgs <- append(lgs, nrow(matched_exclude))
}

pvs_exclude <- data.frame(pvs,lgs)
write.csv(pvs_exclude,"exclusion_at_pvalues.csv", row.names = F)
