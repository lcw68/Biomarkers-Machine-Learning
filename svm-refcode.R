setwd("D:/UNC/research/Xiaojing/antibody")
library(tidyverse)
library(glmnet)
library(caret)
source("svm-ref.R")
X = read.csv("./result/X_ascender_continous.csv",sep=",",header = TRUE)
XX = read.csv("./result/X_reinfect_continous.csv",sep=",",header = TRUE,row.names = 1)
X.binary <- read.csv("./result/X_asc_binary.csv",sep=",",header = TRUE,row.names = 1)
XX.binary <- read.csv("./result/X_reinfect_binary.csv",sep=",",header = TRUE,row.names = 1)

XX.binary <- XX
XX.binary[,2:122] = X.binary[,2:122]
write.csv(XX.binary,file = "./result/X_reinfect_binary.csv")
X[,1] = as.factor(X[,1])
XX[,1] = as.factor(XX[,1])
X.binary[,1] =  as.factor(X.binary[,1])
XX.binary[,1] =  as.factor(XX.binary[,1])

X.binary.filtered <- X.binary[,c(1,which(colSums(X.binary[,-1]) > 8)+1)]
XX.binary.filtered <- XX.binary[,c(1,which(colSums(XX.binary[,-1]) > 8)+1)]

feature_calculate <- function(i,X)
{
  
  #smp_size <- floor(1/3 * nrow(X))
  set.seed(232+7*i)
  train_ind <- createDataPartition(X$yinfect,p=2/3,list = FALSE)[,1]
  nrows = length(train_ind)
  
  nfold = nrows
  folds = rep(1:nfold, len=nrows)[sample(nrows)]
  folds = lapply(1:nfold, function(x) which(folds == x)) ##each training data applied 5-fold svm-rfe
  set.seed(1495)
  results = lapply(folds, svmRFE.wrap, X[train_ind,], k=1, halve.above=30)
  top.features = WriteFeatures(results, X[train_ind,], save=FALSE)
  x = top.features[order(top.features$FeatureID),c("AvgRank")]
  names(x) = colnames(X[,-1])
  return(x)
  #res1[[i]] <- top.features
  #ls <- svmRFE.wrap(test.fold = test_ind, X)
  #rak[,i] = ls$feature.ids
}
#rownames(rak) = colnames(X[,-1])
#head(sort(rowMeans(rak)),10)
svm_cv_error <- function(j,X,select)
{
  set.seed(232+7*j)
  X = na.omit(X)
  train_ind <- createDataPartition(X$yinfect,p=2/3,list = FALSE)[,1]
  #formula0 <- formula(paste("yinfect~",paste(select, collapse=" + ")))
  set.seed(232+7*j)
  tuned = tune.svm(yinfect~.,data = X[train_ind,select],gamma = 10^(-3:-1),cost = 10^(0:3), tunecontrol=tune.control(cross=5))
  return(tuned$best.performance)
}

feature_num <- function(i,j,X,sp)
{
  fetr_select <- c("yinfect",names(head(sp,i)))
  svm_cv_error(j,X,fetr_select)
}

rate_calculate.svm <- function(i,X1,sel)
{
  set.seed(232+7*i)
  train_ind <- createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.use <- X1[train_ind,sel]
  X.test <- X1[-train_ind,sel]
  Y.test <- X.test[,1]
  
  model <- tune.svm(yinfect~.,data = X.use,gamma = 10^(-3:-1),cost = 10^(0:3))
  model0<- svm(yinfect~.,data = X.use,gamma = as.numeric(model$best.parameters[1]),cost = as.numeric(model$best.parameters[2]))
  pr = predict(model0,newdata = X.test)
  #Y.test <- X[test_ind,1]
  # pred = attr(pr,"probabilities")[,2]
  # strc <- roc(Y.test,pred)
  # opt.thre = as.numeric(coords(strc, "best", ret = "threshold",transpose = TRUE))
  # pred.class <- ifelse(pred > opt.thre,1,0)
  dt <- table(factor(pr,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
  T = confusionMatrix(dt)
  return(c(T$overall[1],T$byClass[1:2]))
  
  #gene.error[i] <- mean((X[test_ind,1]-predata)^2)
}

rate1_calculate.svm <- function(i,X1,sel)
{
  set.seed(232+7*i)
  train_ind <- createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.use <- X1[train_ind,sel]
  X.test <- X1[-train_ind,sel]
  Y.test <- X.test[,1]
  
  model <- tune.svm(yinfect~.,data = X.use,gamma = 10^(-3:-1),cost = 10^(0:3))
  model0<- svm(yinfect~.,data = X.use,gamma = as.numeric(model$best.parameters[1]),cost = as.numeric(model$best.parameters[2]),probability = TRUE)
  pr = predict(model0,newdata = X.test,probability=TRUE)
  pred = attr(pr,"probabilities")[,2]
  strc <- roc(Y.test,pred)
  opt.thre = as.numeric(coords(strc, "best", ret = "threshold",transpose = TRUE))
  pred.class <- ifelse(pred > opt.thre,1,0)
  dt <- table(factor(pr,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
  T = confusionMatrix(dt)
  return(c(T$overall[1],T$byClass[1:2]))
  
}

#####main function#####

  rank.ave.asc = sort(rowMeans(do.call(cbind,lapply(1:10,feature_calculate,X.filtered))),decreasing = FALSE)
  rank.ave.asc.bin= sort(rowMeans(do.call(cbind,lapply(1:10,feature_calculate,X.binary.filtered))),decreasing = FALSE)
  rank.ave.reinf= sort(rowMeans(do.call(cbind,lapply(1:10,feature_calculate,XX.filtered%>%dplyr::select(-"prt")))),decreasing = FALSE)
  rank.ave.reinf.bin= sort(rowMeans(do.call(cbind,lapply(1:10,feature_calculate,XX.binary.filtered1))),decreasing = FALSE)
  rank.ave.uninf= sort(rowMeans(do.call(cbind,lapply(1:10,feature_calculate,XXX.filtered))),decreasing = FALSE)
  rank.ave.uninf.bin= sort(rowMeans(do.call(cbind,lapply((1:10)[-7],feature_calculate,XXX.binary.filtered))),decreasing = FALSE)
  
  cv_error.asc <- sapply(1:50,function(i){sapply(1:10,function(j){feature_num(i,j,X.filtered,rank.ave.asc)})})
  cv_error.asc.bin <- sapply(1:50,function(i){sapply(1:10,function(j){feature_num(i,j,X.binary.filtered,rank.ave.asc.bin)})})
  cv_error.reinf <- sapply(1:50,function(i){sapply(1:10,function(j){feature_num(i,j,XX.filtered,rank.ave.reinf)})})
  cv_error.reinf.bin <- sapply(1:50,function(i){sapply(1:10,function(j){feature_num(i,j,XX.binary.filtered,rank.ave.reinf.bin)})})
  cv_error.uninf <- sapply(1:50,function(i){sapply(1:10,function(j){feature_num(i,j,XXX.filtered,rank.ave.uninf)})})
  cv_error.uninf.bin <- sapply(1:44,function(i){sapply(1:10,function(j){feature_num(i,j,XXX.binary.filtered,rank.ave.uninf.bin)})})
  
  rank.ave.reinf.bin = rowMeans(do.call(cbind,lapply(1:s,feature_calculate,XX.binary.filtered)))
  
  
  fetr_select.asc <- c("yinfect",names(head(rank.ave.asc,which.min(colMeans(cv_error.asc)))))
  fetr_select.asc.bin <- c("yinfect",names(head(rank.ave.asc.bin,which.min(colMeans(cv_error.asc.bin)))))
  fetr_select.reinf <- c("yinfect",names(head(rank.ave.reinf,which.min(colMeans(cv_error.reinf)))))
  fetr_select.reinf.bin <- c("yinfect",names(head(rank.ave.reinf.bin,which.min(colMeans(cv_error.reinf.bin)))))
  fetr_select.uninf <- c("yinfect",names(head(rank.ave.uninf,which.min(colMeans(cv_error.uninf)))))
  fetr_select.uninf.bin <- c("yinfect",names(head(rank.ave.uninf.bin,which.min(colMeans(cv_error.uninf.bin)))))
  
  res.asc.con.svm <- do.call(rbind,lapply(1:10,rate_calculate.svm,X.filtered,fetr_select.asc))
  res.asc.bin.svm <- do.call(rbind,lapply(1:10,rate_calculate.svm,X.binary.filtered,fetr_select.asc.bin))
  
  res.reinf.con.svm<- do.call(rbind,lapply(1:10,rate_calculate.svm,XX.filtered,fetr_select.reinf))
  res.reinf.bin.svm<- do.call(rbind,lapply(1:10,rate_calculate.svm,XX.binary.filtered,fetr_select.reinf.bin))
  res.uninf.con.svm<- do.call(rbind,lapply(1:10,rate_calculate.svm,XXX.filtered,fetr_select.uninf))
  res.uninf.bin.svm<- do.call(rbind,lapply(1:10,rate_calculate.svm,XXX.binary.filtered,fetr_select.uninf.bin))
  
