library(caret)
library(tidyverse)
library(ggplot2)
library(pROC)

getwd()
load(file = "antibodydatacode.RData")
setwd("C:/Users/nehav/Local Documents/Research Lab/AntibodyData")
############# DATA ###################
# Ascender Binary
setwd("C:/Users/nehav/Local Documents/Research Lab/AntibodyData")
asc.bin <- read.csv("X_asc_binary.csv",sep=",",header = TRUE)
asc.bin <- na.omit(asc.bin) 
id <- apply(asc.bin[,-1],2,function(x){sum(x!=0,na.rm=TRUE)})
asc.bin.filt <- asc.bin[,c("yinfect",names(id[id>0.05*nrow(asc.bin)]))]  ###keep those number of expression > 5%*total sample
asc.bin.filt[,1]=as.factor(asc.bin.filt[,1])

# Ascender Continous
asc.con <- read.csv("X_ascender_continous.csv",sep=",",header = TRUE)
asc.con.filt <- na.omit(asc.con[,colnames(asc.bin.filt)])###keep the same variables in continuous data
asc.con.filt[,1]=as.factor(asc.con.filt[,1])

# Reinfection Binary
reinf.bin <- read.csv("reinfect_binary_updated.csv",sep=",",header = TRUE)
reinf.bin <- na.omit(reinf.bin) 
id <- apply(reinf.bin[,-1],2,function(x){sum(x!=0,na.rm=TRUE)})
reinf.bin.filt <- reinf.bin[,c("yinfect",names(id[id>0.05*nrow(reinf.bin)]))]  ###keep those number of expression > 5%*total sample
reinf.bin.filt[,1]=as.factor(reinf.bin.filt[,1])

# Reinfection Continous
reinf.con <- read.csv("reinfect_continous_updated.csv",sep=",",header = TRUE)
reinf.con.filt <- na.omit(reinf.con[,colnames(reinf.bin.filt)])###keep the same variables in continuous data
reinf.con.filt[,1]=as.factor(reinf.con.filt[,1])

# Uninfection Binary
uninf.bin <- read.csv("uninfect-binarydata.csv",sep=",",header = TRUE)
uninf.bin <- na.omit(uninf.bin) 
id <- apply(uninf.bin[,-1],2,function(x){sum(x!=0,na.rm=TRUE)})
uninf.bin.filt <- uninf.bin[,c("reinfect",names(id[id>0.05*nrow(uninf.bin)]))]  ###keep those number of expression > 5%*total sample
uninf.bin.filt[,1]=as.factor(uninf.bin.filt[,1])

# Uninfection Continous
uninf.con <- read.csv("uninfect_continuousdata_updated.csv",sep=",",header = TRUE)
uninf.con.filt <- na.omit(uninf.con[,colnames(uninf.bin.filt)])###keep the same variables in continuous data
uninf.con.filt[,1]=as.factor(uninf.con.filt[,1])
uninf.con.filt <- na.omit(uninf.con.filt)

X1 <- uninf.con.filt

knn.varimp <- function(j,X1)   ####calculate variable importance
{
  set.seed(232+7*j)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.train = data.frame(X1[train_ind,-1])
  Y.train = as.factor(X1[train_ind,1])
  trcon = trainControl(method="cv", number=10,search="grid")
  grid <- expand.grid(k = c(5, 11, 21, 25))
  knnmod = train(X.train, Y.train, method= "knn", trControl=trcon,tuneGrid = grid)
  nbs = varImp(knnmod)$importance
  x = nbs[,1]
  names(x) = rownames(nbs)
  return(x)
}

knn.var.rank.asc.con = sort(rowMeans(do.call(cbind,lapply(1:10,knn.varimp,X.filtered))),decreasing = TRUE)
knn.var.rank.asc.bin= sort(rowMeans(do.call(cbind,lapply(1:10,knn.varimp,X.binary.filtered))),decreasing = TRUE)

knn.var.rank.reinf.con = sort(rowMeans(do.call(cbind,lapply(1:10,knn.varimp,XX.filtered))),decreasing = TRUE)
knn.var.rank.reinf.bin = sort(rowMeans(do.call(cbind,lapply(1:10,knn.varimp,XX.binary.filtered))),decreasing = TRUE)

knn.var.rank.uninf.bin = sort(rowMeans(do.call(cbind,lapply(1:10,knn.varimp,XXX.binary.filtered))),decreasing = TRUE)
knn.var.rank.uninf.con = sort(rowMeans(do.call(cbind,lapply(1:10,knn.varimp,XXX.filtered))),decreasing = TRUE)


num.calculate.knn <- function(i,sel,X1) ### calculate best number in the final model
{
  select.var <- names(sel)[1:i]
  accu = 0
  for(j in 1:10)
  {
    set.seed(232+7*j)
    train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
    X.train.use <- X1[train_ind,select.var]
    Y.train.use = as.factor(X1[train_ind,1])
    X.use <- X1[train_ind,c("yinfect",select.var)]
    trcon = trainControl(method="cv", number=10,search="grid",classProbs = TRUE,summaryFunction = twoClassSummary)
    grid <- expand.grid(k = c(5, 11, 21, 25))
    X.use <- X.use %>% mutate(yinfect = factor(yinfect, 
                                               labels = make.names(levels(yinfect))))
    knnmod = train(yinfect~.,data=X.use,method="knn",trControl = trcon,tuneGrid = grid,metric='ROC')
    accu[j] <- max(knnmod$results[,"ROC"],na.rm=TRUE)
  }
  perfom <- mean(accu,na.rm=TRUE)
  
}

knn.perform.res.asc.con = unlist(lapply(2:35,num.calculate.knn,knn.var.rank.asc.con,X.filtered))
knn.perform.res.asc.bin = unlist(lapply(2:35,num.calculate.knn,knn.var.rank.asc.bin,X.binary.filtered))

knn.perform.res.reinf.con = unlist(lapply(2:35,num.calculate.knn,knn.var.rank.reinf.con,XX.filtered))
knn.perform.res.reinf.bin = unlist(lapply(2:35,num.calculate.knn,knn.var.rank.reinf.bin,XX.binary.filtered))
knn.perform.res.uninf.con = unlist(lapply(2:35,num.calculate.knn,knn.var.rank.uninf.con,XXX.filtered))

knn.perform.res.uninf.bin = unlist(lapply(2:35,num.calculate.knn,knn.var.rank.uninf.bin,XXX.binary.filtered))


final.select.knn.asc.bin <- names(knn.var.rank.asc.bin)[1:which.max(knn.perform.res.asc.bin)]

final.select.knn.asc.con <- names(knn.var.rank.asc.con)[1:which.max(knn.perform.res.asc.con)]

final.select.knn.reinf.bin <- names(knn.var.rank.reinf.bin)[1:which.max(knn.perform.res.reinf.bin)]

final.select.knn.reinf.con <- names(knn.var.rank.reinf.con)[1:which.max(knn.perform.res.reinf.con)]

final.select.knn.uninf.bin <- names(knn.var.rank.uninf.bin)[1:which.max(knn.perform.res.uninf.bin)]

final.select.knn.uninf.con <- names(knn.var.rank.uninf.con)[1:which.max(knn.perform.res.uninf.con)]


rate_calculate.knn <- function(j,select.var,X1) 
{
  set.seed(232+7*j)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.test <- X1[-train_ind,select.var]
  Y.test = as.factor(X1[-train_ind,1])
  
  X.use <- X1[train_ind,c("yinfect",select.var)]
  X.use <- X.use %>% mutate(yinfect = factor(yinfect,labels = make.names(levels(yinfect))))
  grid <- expand.grid(k = c(5, 11, 21, 25))
  trcon = trainControl(classProbs = TRUE,method="cv",number = 10, search="grid",summaryFunction = twoClassSummary)
  trmod = train(yinfect~.,data=X.use,method="knn",trControl = trcon,metric='ROC',tuneGrid = grid)
  pred = predict(trmod,newdata = X.test,type="prob")
  rocm <- pROC::roc(Y.test,pred[,2])
  opt.thre = as.numeric(coords(rocm, "best", ret = "threshold",transpose = TRUE))
  pred.class <- ifelse(pred[,2] > opt.thre,1,0)
  dt <- table(factor(pred.class,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
  T = confusionMatrix(dt)
  return(c(T$overall[1],T$byClass[1:2]))
}


res.asc.con.knn = do.call(rbind,lapply(1:10,rate_calculate.knn,final.select.knn.asc.con,X.filtered))
res.asc.bin.knn = do.call(rbind,lapply(1:10,rate_calculate.knn,final.select.knn.asc.bin,X.binary.filtered))

res.reinf.con.knn = do.call(rbind,lapply(1:10,rate_calculate.knn,final.select.knn.reinf.con,XX.filtered))
res.reinf.bin.knn = do.call(rbind,lapply(1:10,rate_calculate.knn,final.select.knn.reinf.bin,XX.binary.filtered))

res.uninf.con.knn= do.call(rbind,lapply(1:10,rate_calculate.knn,final.select.knn.uninf.con,XXX.filtered))

res.uninf.bin.knn = do.call(rbind,lapply(1:10,rate_calculate.knn,final.select.knn.uninf.bin,XXX.binary.filtered))








