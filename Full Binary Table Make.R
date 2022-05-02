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

#set up matched snps
matched_snps_file <- read.table("matched_snps_annotated.txt", sep = '\t', header = T)
matched_snps <- matched_snps_file %>% select('snpID', 'rsID', 'snp_maf', 
                                             'dist_nearest_gene', 'friends_ld07')
excluded_matches <- read.csv('excluded_matches.csv', header = T)
matched_snps <- matched_snps[!(matched_snps$rsID %in% excluded_matches$rsID),]
matched_snps <- matched_snps %>% separate(snpID, c("chr", "pos"), remove = F)
matched_snps$chr <- paste0('chr', matched_snps$chr)
matched_snps$GWAS <- 0

#set up input snps
input_snps_file <- read.table("input_snps_annotated.txt", sep = '\t', header = T)
input_snps <- input_snps_file %>% select('snpID', 'rsID', 'snp_maf', 
                                         'dist_nearest_gene', 'friends_ld07')
input_snps <- input_snps %>% separate(snpID, c("chr", "pos"), remove = F)
input_snps$chr <- paste0('chr', input_snps$chr)
input_snps$GWAS <- 1

#write bed file of combined snps
bed_file <- rbind(matched_snps %>% select('chr','pos'), input_snps %>% select('chr', 'pos'))
write.table(bed_file, file = "snps_matched_.bed", 
            row.names=F, sep="\t", col.names = F)

#create master binary table
binary_table <- rbind(input_snps, matched_snps)

#load file with names of factors and files
files <- list.files(path="/Users/kamran/Documents/Spring 2022/Research/SNP GWAS Data/Sample Encode Data/", pattern="*.", full.names=TRUE, recursive=FALSE)

for (i in files) {
  factor_list <- c()
  factor_bed_ <- fread(i, select = c(1:3))
  colnames(factor_bed_) <- c('chr', 'start', 'stop')
  name <- gsub("\\..*","",basename(i))
  print(name)
  for (j in 1:nrow(binary_table)) {
    factor_bed <- factor_bed_
    factor_bed <- factor_bed[strtoi(factor_bed$start) <= strtoi(binary_table[j,'pos']),]
    factor_bed <- factor_bed[strtoi(factor_bed$stop) >= strtoi(binary_table[j,'pos']),]
    factor_bed <- factor_bed[factor_bed$chr == binary_table[j,'chr'],]
    if (nrow(factor_bed) > 0) {
      factor_list <- append(factor_list, 1)
    }
    else {
      factor_list <- append(factor_list, 0)
    }
  }
  binary_table <- cbind(binary_table, factor_list)
  colnames(binary_table)[ncol(binary_table)] <- name
}

write.csv(binary_table, 'binary_table.csv', row.names = F)