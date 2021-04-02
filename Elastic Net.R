####

library(tidyverse)
library(glmnet)
library(caret)
source("svm-ref.R")
X = read.csv("X_ascender_continous.csv",sep=",",header = TRUE,row.names = 1)
XX = read.csv("X_reinfect_continous.csv",sep=",",header = TRUE,row.names = 1)
X.binary <- read.csv("X_asc_binary.csv",sep=",",header = TRUE,row.names = 1)
XX.binary <- XX
XX.binary[,2:122] = X.binary[,2:122]

X[,1] = as.factor(X[,1])
XX[,1] = as.factor(XX[,1])
X.binary[,1] =  as.factor(X.binary[,1])
X.binary.filtered <- X.binary[,c(1,which(colSums(X.binary[,-1]) > 8)+1)]
XX.binary.filtered <- XX.binary[,c(1,which(colSums(XX.binary[,-1]) > 8)+1)]

get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

els.varimp <- function(i,X1)
{
  set.seed(232+7*i)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.train = as.matrix(X1[train_ind,-1])
  Y.train = X1[train_ind,1]
  trctrl = trainControl(method = "cv", number = 10, search ="grid")
  alpha  <- seq(0,1,length.out = 25)
  lambda <- seq(0,1,length.out = 25)
  #set.seed(1234)
  #for(k in 1:10)
  #{
  #set.seed(2345+k)
  set.seed(232+7*i)
  
  dw_elnet <- train(X.train, Y.train, 
                    method = "glmnet", trControl = trctrl, family = "binomial", 
                    tuneGrid = expand.grid(alpha = alpha, lambda = lambda))
  #rept.cv <- rbind(rept.cv,get_best_result(dw_elnet))
  #}
  imp.mat <- varImp(dw_elnet,lambda = dw_elnet$bestTune$lambda,alpha = dw_elnet$bestTune$alpha)$importance
  x = imp.mat[,1]
  names(x) = rownames(imp.mat)
  return(x)
}

els.var.rank.asc.con = sort(rowMeans(do.call(cbind,lapply(1:10,els.varimp,X.filtered)),na.rm=TRUE),decreasing = TRUE)
els.var.rank.asc.bin = sort(rowMeans(do.call(cbind,lapply(1:10,els.varimp,X.binary.filtered)),na.rm=TRUE),decreasing = TRUE)
els.var.rank.reinf.con = sort(rowMeans(do.call(cbind,lapply(1:10,els.varimp,XX.filtered)),na.rm=TRUE),decreasing = TRUE)
els.var.rank.reinf.bin = sort(rowMeans(do.call(cbind,lapply(1:10,els.varimp,XX.binary.filtered)),na.rm=TRUE),decreasing = TRUE)
els.var.rank.uninf.con = sort(rowMeans(do.call(cbind,lapply(1:10,els.varimp,XXX.filtered)),na.rm=TRUE),decreasing = TRUE)
els.var.rank.uninf.bin = sort(rowMeans(do.call(cbind,lapply(1:10,els.varimp,XXX.binary.filtered)),na.rm=TRUE),decreasing = TRUE)

elas_sp_asc <- sort(table(unlist(elas_sel_asc)),decreasing =TRUE)

model0 <- glmnet(X.train, Y.train, family = "binomial", 
                 alpha = dw_elnet$bestTune$alpha, lambda = dw_elnet$bestTune$lambda)
coefm <- as.matrix((coef(model0)))
varselected <- names(coefm[which(coefm[,1] !=0),1])[-1]
num1 <- sum(elas_sp_asc > 5)
X1<- na.omit(X)

num.calculate.els<- function(i,sel,X1)
{
  spp <- names(sel[1:i])
  performance = 0
  for(j in 1:10)
  {
    set.seed(232+7*j)
    train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
    # X.train = as.matrix(X1[train_ind,spp])
    # Y.train = X1[train_ind,1]
    X.use <- X1[train_ind,c("yinfect",spp)]
    cv_5 = trainControl(method = "cv", number = 5,search = "grid")
    alpha  <- seq(0,1,length.out = 25)
    lambda <- seq(0,1,length.out = 25)
    set.seed(232+7*j)
    # X.use <- X.use %>% mutate(yinfect = factor(yinfect, 
    #                                            labels = make.names(levels(yinfect))))
    dw_elnet <- train(yinfect~.,data=X.use, 
                      method = "glmnet", trControl = cv_5, family = "binomial", 
                      tuneGrid = expand.grid(alpha = alpha, lambda = lambda))
    performance[j] <- get_best_result(dw_elnet)[3]
  }
  perform = mean(unlist(performance),na.rm=TRUE)
  return(perform)
  
}
els.perform.res.asc.con = unlist(lapply(2:30,num.calculate.els,els.var.rank.asc.con,X.filtered))

els.perform.res.asc.bin = unlist(lapply(2:30,num.calculate.els,els.var.rank.asc.bin,X.binary.filtered))
els.perform.res.reinf.con = unlist(lapply(2:32,num.calculate.els,els.var.rank.reinf.con,XX.filtered))

els.perform.res.reinf.bin = unlist(lapply(2:32,num.calculate.els,els.var.rank.reinf.bin,XX.binary.filtered))

els.perform.res.uninf.con = unlist(lapply(2:35,num.calculate.els,els.var.rank.uninf.con,XXX.filtered))

els.perform.res.uninf.bin = unlist(lapply(2:35,num.calculate.els,els.var.rank.uninf.bin,XXX.binary.filtered))
els.final.select.asc.con <- names(els.var.rank.asc.con)[1:(which.max(els.perform.res.asc.con)+1)]
els.final.select.asc.bin <- names(els.var.rank.asc.bin)[1:(which.max(els.perform.res.asc.bin)+1)]

els.final.select.reinf.con <- names(els.var.rank.reinf.con)[1:(which.max(els.perform.res.reinf.con)+1)]
els.final.select.reinf.bin <- names(els.var.rank.reinf.bin)[1:(which.max(els.perform.res.reinf.bin)+1)]
els.final.select.uninf.bin <- names(els.var.rank.uninf.bin)[1:(which.max(els.perform.res.uninf.bin)+1)]
els.final.select.uninf.con <- names(els.var.rank.uninf.con)[1:(which.max(els.perform.res.uninf.con)+1)]


spp.asc <- names(elas_sp_asc[1:which.max(acc.asc)])

rate_calculate.els<- function(j,spp,X1)
{
  set.seed(132+7*j)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.train = as.matrix(X1[train_ind,spp])
  Y.train = X1[train_ind,1]
  X.test <- as.matrix(X1[-train_ind,spp])
  Y.test = X1[-train_ind,1]
  cv_5 = trainControl(method = "cv", number = 10)
  alpha  <- seq(0,1,length.out = 50)
  lambda <- seq(0,1,length.out = 50)
  set.seed(232+34*j)
  dw_elnet <- train(X.train, Y.train, 
                    method = "glmnet", trControl = cv_5, family = "binomial", 
                    tuneGrid = expand.grid(alpha = alpha, lambda = lambda))
 
  model0 <- glmnet(X.train, Y.train, family = "binomial", 
                   alpha = dw_elnet$bestTune$alpha, lambda = dw_elnet$bestTune$lambda)

  pred <- predict(model0, newx = X.test, type = "response",
                     alpha = dw_elnet$bestTune$alpha, lambda = dw_elnet$bestTune$lambda)
  rocm <- roc(Y.test,pred[,1])
  opt.thre = as.numeric(coords(rocm, "best", ret = "threshold",transpose = TRUE))
  pred.class <- ifelse(pred[,1] > opt.thre,1,0)
  dt <- table(factor(pred.class,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
  T = confusionMatrix(dt)
  return(c(T$overall[1],T$byClass[1:2]))
  #predata = ifelse(pred >0.5,1,0)
  
}

res.asc.con.els = do.call(rbind,lapply(1:10,rate_calculate.els,els.final.select.asc.con,X.filtered))
res.asc.bin.els = do.call(rbind,lapply(1:10,rate_calculate.els,els.final.select.asc.bin,X.binary.filtered))

res.reinf.con.els = do.call(rbind,lapply(1:10,rate_calculate.els,els.final.select.reinf.con,XX.filtered))
res.reinf.bin.els = do.call(rbind,lapply(1:10,rate_calculate.els,els.final.select.reinf.bin,XX.binary.filtered))

res.uninf.con.els1= do.call(rbind,lapply(1:10,rate_calculate.els,els.final.select.uninf.con,XXX.filtered))

res.uninf.bin.els = do.call(rbind,lapply(1:10,rate_calculate.els,els.final.select.uninf.bin,XXX.binary.filtered))


els.final.select.asc.bin.inter = intersect(els.final.select.asc.bin,protein)
els.final.select.reinf.bin.inter = intersect(els.final.select.reinf.bin,protein)
els.final.select.uninf.bin.inter = intersect(els.final.select.uninf.bin,protein)

rocplot <- function(j,spp,X1)
{
  set.seed(132+7*j)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.train = as.matrix(X1[train_ind,spp])
  Y.train = X1[train_ind,1]
  X.test <- as.matrix(X1[-train_ind,spp])
  Y.test = X1[-train_ind,1]
  cv_5 = trainControl(method = "cv", number = 10)
  alpha  <- seq(0,1,length.out = 50)
  lambda <- seq(0,1,length.out = 50)
  set.seed(232+34*j)
  dw_elnet <- train(X.train, Y.train, 
                    method = "glmnet", trControl = cv_5, family = "binomial", 
                    tuneGrid = expand.grid(alpha = alpha, lambda = lambda))
  model0 <- glmnet(X.train, Y.train, family = "binomial", 
                   alpha = dw_elnet$bestTune$alpha, lambda = dw_elnet$bestTune$lambda)
  
  pred <- predict(model0, newx = X.test, type = "response",
                  alpha = dw_elnet$bestTune$alpha, lambda = dw_elnet$bestTune$lambda)
  modr = roc(Y.test,as.numeric(pred))
  return(modr)
}

plot(rocplot(5,els.final.select.asc.con,X.filtered),print.auc=TRUE,auc.polygon=TRUE,grid=c(0.1,0.2),grid.col=c("green","red"),max.auc.polygon=TRUE,auc.polygon.col="skyblue",#print.thres=TRUE,
     main = "")

