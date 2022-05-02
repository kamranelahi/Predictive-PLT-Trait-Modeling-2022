binary_table <- read.csv('binary_table_wblood.csv', header = T)
setwd('./Lasso Plots/With Interactions/')
cor_binary <- cor(binary_table[,-c(1:8)], method = c("pearson", "kendall", "spearman"))

test_limit <- function(cutoff){
  n <- c()
  for (i in 1:ncol(cor_binary)){
    for (j in 1:ncol(cor_binary)){
      if (cor_binary[i,j] > cutoff && j!=i){
        n <- append(n, colnames(cor_binary)[i])
      }
    }
  }
  z <- binary_table[,!names(binary_table) %in% n]
  
  x <- data.matrix(z[,c(5:7,9:ncol(z))])
  y <- data.matrix(z[,8])
  n <- ncol(x)
  pfactor = matrix(nrow = 1, ncol = n)
  pfactor[1,] = 1
  pfactor[1,c(1:3)] = 0
  cvfit = cv.glmnet(x, y, penalty.factor = pfactor, family = "binomial", type.measure = "auc")
  
  pdf(paste0("subset",toString(cutoff), "LASSO.output.AUCplot.pdf"))
  plot(cvfit, main = toString(cutoff))
  dev.off()
  
  cat("coef(cvfit, s=lambda.1se)\n")
  coef(cvfit, s="lambda.1se")
  cat("coef(cvfit, s=lambda.min)\n")
  coef(cvfit, s="lambda.min")
}

for (i in seq(.5,1, .1)){
  test_limit(i)
  print(i)
}
