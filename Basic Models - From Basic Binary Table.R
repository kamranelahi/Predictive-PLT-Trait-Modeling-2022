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

#prepare table for analysis
gene_factors <- binary_table %>% select(-c('snpID','chr','pos'))
shuffled_data= gene_factors[sample(1:nrow(gene_factors)), ]
rsIDs <- shuffled_data$rsID
shuffled_data$rsID <- NULL
cutoff <- as.integer(0.8*nrow(shuffled_data))
train_factors <- shuffled_data[1:cutoff,]
test_factors <- shuffled_data[(cutoff+1):nrow(shuffled_data),]
GWAS_test <- test_factors$GWAS
test_factors$GWAS <- NULL

#simple logistic regression
glm_fit = glm(GWAS ~ (.),  family = "binomial", data = train_factors)
summary(glm_fit)
fitted_probabilities_glm = predict(glm_fit, newdata = test_factors, type = "response")

auc(roc(GWAS_test,fitted_probabilities_glm))

#lasso penalized regression
x_train_matrix = model.matrix( GWAS ~ ., data = train_factors)[,-1]  
y_train        = train_factors$GWAS

x_test_matrix  = model.matrix( GWAS ~ ., data = test_factors )[,-1]  
y_test         = test_factors$GWAS

lasso_model = cv.glmnet(x = x_train_matrix, y = as.factor(y_train), alpha = 1, family = 'binomial') 
lasso_model
lasso_predictions = predict(lasso_model, newx = x_test_matrix, s = lasso_model$lambda.1se)

auc(roc(GWAS_test,lasso_predictions))

#random forest
shuffled_data_ = shuffled_data %>%
  mutate(GWAS = factor(GWAS)) 
n_total = nrow(shuffled_data)
n_train = round(0.8*n_total)
n_test = n_total-n_train
partition = sample(c(rep("train", n_train), rep("test", n_test)))
factor_train = shuffled_data_ %>%
  bind_cols(partition = partition) %>%
  filter(partition == "train") %>%
  select(-partition)
factor_test = shuffled_data_ %>%
  bind_cols(partition = partition) %>%
  filter(partition == "test") %>%
  select(-partition)

names(factor_train) <- gsub(" ", "_", names(factor_train))
names(factor_train) <- gsub("[()]", "_", names(factor_train))
names(factor_train) <- gsub("[+]", "", names(factor_train))
names(factor_train) <- gsub("[-]", "", names(factor_train))

names(factor_test) <- gsub(" ", "_", names(factor_train))
names(factor_test) <- gsub("[()]", "_", names(factor_train))
names(factor_test) <- gsub("[+]", "", names(factor_train))
names(factor_test) <- gsub("[-]", "", names(factor_train))

rf_fit = randomForest(GWAS ~ ., data = factor_train)
# rf_fit$importance
# varImpPlot(rf_fit)

pdf('rf_plot.pdf')
tibble(oob_error = rf_fit$err.rate[,"OOB"],
       trees = 1:500) %>%
  ggplot(aes(x = trees, y = oob_error)) + geom_line() + theme_bw()
dev.off()

rf_predictions = predict(rf_fit, newdata = factor_test)

tmp_rf <- as.numeric(levels(rf_predictions))[rf_predictions]

auc(roc(GWAS_test,tmp_rf))

#boosting

#train_factors$GWAS <- as.numeric(levels(train_factors$GWAS))[train_factors$GWAS]
gbm_fit_1 = gbm(GWAS ~ .,
                distribution = "bernoulli",
                n.trees = 2000,
                interaction.depth = 1,
                shrinkage = 0.01,
                cv.folds = 5,
                data = train_factors)
gbm_fit_2 = gbm(GWAS ~ .,
                distribution = "bernoulli",
                n.trees = 2000,
                interaction.depth = 2,
                shrinkage = 0.01,
                cv.folds = 5,
                data = train_factors)
gbm_fit_3 = gbm(GWAS ~ .,
                distribution = "bernoulli",
                n.trees = 2000,
                interaction.depth = 3,
                shrinkage = 0.01,
                cv.folds = 5,
                data = train_factors)
gbm_fit_4 = gbm(GWAS ~ .,
                distribution = "bernoulli",
                n.trees = 2000,
                interaction.depth = 4,
                shrinkage = 0.01,
                cv.folds = 5,
                data = train_factors)
gbm_fit_5 = gbm(GWAS ~ .,
                distribution = "bernoulli",
                n.trees = 2000,
                interaction.depth = 5,
                shrinkage = 0.01,
                cv.folds = 5,
                data = train_factors)

ntrees = 2000
cv_errors = bind_rows(
  tibble(ntree = 1:ntrees, cv_err = gbm_fit_1$cv.error, depth = 1),
  tibble(ntree = 1:ntrees, cv_err = gbm_fit_2$cv.error, depth = 2),
  tibble(ntree = 1:ntrees, cv_err = gbm_fit_3$cv.error, depth = 3),
  tibble(ntree = 1:ntrees, cv_err = gbm_fit_4$cv.error, depth = 4),
  tibble(ntree = 1:ntrees, cv_err = gbm_fit_5$cv.error, depth = 5)
)

pdf('boosting_errors_12345.pdf')
cv_errors %>%
  ggplot(aes(x = ntree, y = cv_err, colour = factor(depth))) +
  geom_line() + theme_bw()
dev.off()

gbm_fit_optimal = gbm_fit_2
optimal_num_trees = gbm.perf(gbm_fit_2, plot.it = FALSE)
summary(gbm_fit_optimal, n.trees = optimal_num_trees, plotit = FALSE)
write.csv(summary(gbm_fit_optimal, n.trees = optimal_num_trees, plotit = FALSE), "boost_perf.csv", row.names = F)

gbm_probabilities = predict(gbm_fit_optimal, n.trees = optimal_num_trees,
                            type = "response", newdata = test_factors)
# gbm_predictions = as.numeric(gbm_probabilities > mean(gbm_probabilities+6*sd(gbm_probabilities)))
# result_comparison <- data.frame(y_test, gbm_predictions)
# result_comparison$correct <- ifelse(result_comparison$y_test == result_comparison$gbm_predictions, 1, 0)
# top <- mean(result_comparison$correct)
# 
# for (i in 1:100){
#   gbm_predictions = as.numeric(gbm_probabilities > mean(gbm_probabilities+i*sd(gbm_probabilities)))
#   result_comparison <- data.frame(y_test, gbm_predictions)
#   result_comparison$correct <- ifelse(result_comparison$y_test == result_comparison$gbm_predictions, 1, 0)
#   if( mean(result_comparison$correct) > top){
#     print(paste(i, mean(result_comparison$correct), sep = ' '))
#   }
# }

auc(roc(GWAS_test,gbm_probabilities))

#averaging factors
lasso_logistic <- c()
for (i in 1:length(fitted_probabilities_glm)) {
  lasso_logistic <- append(lasso_logistic, mean(fitted_probabilities_glm[i], lasso_predictions[i]))
}
auc(roc(GWAS_test,lasso_logistic))

boosting_logistic <- c()
for (i in 1:length(fitted_probabilities_glm)) {
  boosting_logistic <- append(boosting_logistic, mean(fitted_probabilities_glm[i], gbm_probabilities[i]))
}
auc(roc(GWAS_test,boosting_logistic))

boosting_lasso <- c()
for (i in 1:length(gbm_probabilities)) {
  boosting_lasso <- append(boosting_lasso, mean(lasso_predictions[i], gbm_probabilities[i]))
}
auc(roc(GWAS_test,boosting_lasso))

all_avg <- c()
for (i in 1:length(gbm_probabilities)) {
  all_avg <- append(all_avg, mean(lasso_predictions[i], gbm_probabilities[i], fitted_probabilities_glm[i]))
}
auc(roc(GWAS_test,all_avg))
