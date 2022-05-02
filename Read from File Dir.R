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