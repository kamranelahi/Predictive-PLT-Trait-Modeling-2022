z <- read.csv('binary_table_wblood.csv', header=T)
library(glmnet)
sink("GLMnet_Lasso_output.txt") #Note that column 259 was corrupted, hence left out of analysis
cat("\nnrow: ")
nrow(z)
cat("ncol: ")
ncol(z)
x <- data.matrix(z[,c(5:7,9:37)])
y <- data.matrix(z[,8])
pfactor = matrix(nrow = 1, ncol = 32)
pfactor[1,] = 1
pfactor[1,c(1:3)] = 0
cvfit = cv.glmnet(x, y, penalty.factor = pfactor, family = "binomial", type.measure = "auc")

cat("print(cvfit)\n")
print(cvfit)
cat("coef(cvfit, s=lambda.1se)\n")
coef(cvfit, s="lambda.1se")
cat("coef(cvfit, s=lambda.min)\n")
coef(cvfit, s="lambda.min")
cat("AUC saved LASSO.output.AUCplot.pdf\n")
pdf("subsetLASSO.output.AUCplot.pdf")
plot(cvfit)
dev.off()