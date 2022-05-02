binary_table <- read.csv('binary_table_wblood.csv', header = T)
setwd('./Lasso Plots/With Interactions/')
cor_binary <- cor(binary_table[,c(237:254)], method = c("pearson", "kendall", "spearman"))

test_limit <- function(cutoff){
  n <- data.frame(binary_table[,1])
  for (i in 1:ncol(cor_binary)){
    for (j in 1:ncol(cor_binary)){
      if (cor_binary[i,j] > cutoff && j!=i){
        n <- cbind(n, binary_table[,colnames(cor_binary)[i]]*binary_table[,colnames(cor_binary)[j]])
        colnames(n)[ncol(n)] <- paste0(colnames(cor_binary)[i],"*",colnames(cor_binary)[j])
      }
    }
  }
  n <- n[,-1]
  z <- cbind(binary_table,n)
  
  x_ <- data.matrix(z[,c(5:7,9:ncol(z))])
  y <- data.matrix(z[,8])
  n <- ncol(x_)
  pfactor = matrix(nrow = 1, ncol = n)
  pfactor[1,] = 1
  pfactor[1,c(1:3)] = 0
  cvfit = cv.glmnet(x_, y, penalty.factor = pfactor, family = "binomial", type.measure = "auc")
  
  #print(coef(cvfit,s='lambda.1se'))
  
  pdf(paste0("subset",toString(cutoff), "LASSO.output.AUCplot.pdf"))
  plot(cvfit, main = toString(cutoff))
  dev.off()

  cat("coef(cvfit, s=lambda.1se)\n")
  coef(cvfit, s="lambda.1se")
  cat("coef(cvfit, s=lambda.min)\n")
  coef(cvfit, s="lambda.min")
}

for (i in seq(.1, 1, .1)){
  test_limit(i)
}

