library(caret)
library(tidyverse)
library(randomForest)
X1 = na.omit(X.binary)
X1 = na.omit(XX)
X1[,1]=as.factor(X1[,1])
X[,1] = as.factor(X[,1])
XX[,1] = as.factor(XX[,1])
X.binary[,1] =  as.factor(X.binary[,1])
XX.binary[,1] =  as.factor(XX.binary[,1])

X.binary.filtered <- na.omit(X.binary[,c(1,which(colSums(X.binary[,-1],na.rm=TRUE) >= 8)+1)])

X.binary1 <- na.omit(X.binary)
id <- apply(X.binary1[,-1],2,function(x){sum(x!=0,na.rm=TRUE)})
X.binary.filtered <- X.binary1[,c("yinfect",names(id[id>0.05*nrow(X.binary1)]))]
X.filtered <- na.omit(X[,colnames(X.binary.filtered)])

XX.binary1 <- na.omit(XX.binary)
id <- apply(XX.binary1[,-1],2,function(x){sum(x!=0,na.rm=TRUE)})
XX.binary.filtered <- XX.binary1[,c("yinfect",names(id[id>0.05*nrow(XX.binary1)]))]
XX.filtered <- na.omit(XX[,colnames(XX.binary.filtered)])

XXX.binary1 <- na.omit(bin.uninfect)
id <- apply(XXX.binary1[,-1],2,function(x){sum(x!=0,na.rm=TRUE)})
XXX.binary.filtered <- XXX.binary1[,c("yinfect",names(id[id>0.05*nrow(XXX.binary1)]))]
XXX.filtered <- na.omit(con.uninfect[,colnames(XXX.binary.filtered)])

XXX.filtered[,1] = as.factor(XXX.filtered[,1])
XXX.binary.filtered[,1] = as.factor(XXX.binary.filtered[,1])
XX.binary1 <- na.omit(XX.binary)
id <- apply(XX.binary1[,-1],2,function(x){sum(x!=0,na.rm=TRUE)})
XX.binary.filtered <- XX.binary1[,c("yinfect",names(id[id>5]))]

X1 = XXX.binary.filtered
for(j in 1:10)
{
  set.seed(232+7*j)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.train = data.frame(X1[train_ind,-1])
  Y.train = as.factor(X1[train_ind,1])
  X.test = data.frame(X1[-train_ind,-1])
  Y.test = as.factor(X1[-train_ind,1])
  #rfImp1 <- randomForest(yinfect ~ ., data = X1[train_ind,], ntree = 1000, importance = TRUE)
  if(method == "nb")
    
  {
    ctrl <- rfeControl(functions = nbFuncs,
                       method = "cv",number = 10,
                       verbose = FALSE,rerank = FALSE)
  }
  
  else if (method == "rf")
  {
    ctrl <- rfeControl(functions = rfFuncs,
                       method = "cv",number = 10,
                       verbose = FALSE,rerank = FALSE)
  }
    
  rfprofile = rfe(X.train,Y.train,sizes = 1:45,rfeControl = ctrl)
  ls.nb.con[[j]] <- rfprofile$optVariables
  #return(importance(rfImp1)[,3])
  pred1 <- predict(rfprofile,newdata = X.test)$pred
  dt <- table(factor(pred1,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
  T1 = confusionMatrix(dt)
  rfe.res.con.nb <- rbind(rfe.res.con.nb,c(T1$overall[1],T1$byClass[1:2]))
}

Imp_calculate <- function(j,X1)
{
  set.seed(232+7*j)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  rfImp1 <- randomForest(yinfect ~ ., data = X1[train_ind,], ntree = 1000, importance = TRUE)
  return(importance(rfImp1)[,3])
}


var.rank.asc.con = sort(rowMeans(do.call(cbind,lapply(1:10,Imp_calculate,X.filtered))),decreasing = TRUE)
var.rank.asc.bin = sort(rowMeans(do.call(cbind,lapply(1:10,Imp_calculate,X.binary.filtered))),decreasing = TRUE)
var.rank.reinf.con = sort(rowMeans(do.call(cbind,lapply(1:10,Imp_calculate,XX.filtered))),decreasing = TRUE)
var.rank.reinf.bin = sort(rowMeans(do.call(cbind,lapply(1:10,Imp_calculate,XX.binary.filtered))),decreasing = TRUE)
var.rank.uninf.con = sort(rowMeans(do.call(cbind,lapply(1:10,Imp_calculate,XXX.filtered))),decreasing = TRUE)
var.rank.uninf.bin = sort(rowMeans(do.call(cbind,lapply(1:10,Imp_calculate,XXX.binary.filtered))),decreasing = TRUE)


num.calculate <- function(i,X1,var.rank)
{
  select.var <- names(var.rank)[1:i]
  oob.error = 0
  for(j in 1:10)
  {
    set.seed(232+7*j)
    train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
    X.train.use <- X1[train_ind,select.var]
    Y.train.use = as.factor(X1[train_ind,1])
    X.use <- X1[train_ind,c("yinfect",select.var)]
    trmod1 <- tuneRF(X.train.use,Y.train.use,stepFactor = 2,ntreeTry = 50)
    oob.error[j] <- min(trmod1[,2])
  }
  perfom <- mean(oob.error,na.rm=TRUE)
  #tunegrid <- expand.grid(.mtry=1:i)
  #trcon = trainControl(method="repeatedcv", number=10, repeats=3,search="grid")
  #trmod = train(yinfect~.,data=X.use,ntree=1000,tuneGrid = tunegrid,trControl = trcon,metric='Accuracy')
}

perform.res.asc.con = unlist(lapply(2:40,num.calculate,X.filtered,var.rank.asc.con))
perform.res.asc.bin <-  unlist(lapply(2:40,num.calculate,X.binary.filtered,var.rank.asc.bin))
perform.res.reinf.con = unlist(lapply(2:40,num.calculate,XX.filtered,var.rank.reinf.con))
perform.res.reinf.bin <-  unlist(lapply(2:40,num.calculate,XX.binary.filtered,var.rank.reinf.bin))
perform.res.uninf.con = unlist(lapply(2:40,num.calculate,XXX.filtered,var.rank.uninf.con))
perform.res.uninf.bin <-  unlist(lapply(2:40,num.calculate,XXX.binary.filtered,var.rank.uninf.bin))

final.select.asc.con <- names(var.rank.asc.con)[1:(which.min(perform.res.asc.con)+1)]
final.select.asc.bin <- names(var.rank.asc.bin)[1:(which.min(perform.res.asc.bin)+1)]
final.select.reinf.con <- names(var.rank.reinf.con)[1:(which.min(perform.res.reinf.con)+1)]
final.select.reinf.bin <- names(var.rank.reinf.bin)[1:(which.min(perform.res.reinf.bin)+1)]
final.select.uninf.con <- names(var.rank.uninf.con)[1:(which.min(perform.res.uninf.con)+1)]
final.select.uninf.bin <- names(var.rank.uninf.bin)[1:(which.min(perform.res.uninf.bin)+1)]

rate_calculate.rf <- function(j,X1,select.var)
{
  
  set.seed(232+7*j)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.test <- X1[-train_ind,select.var]
  Y.test = as.factor(X1[-train_ind,1])
  X.use <- X1[train_ind,c("yinfect",select.var)]
  tunegrid <- expand.grid(.mtry=seq(1,length(select.var),2))
  trcon = trainControl(method="cv", number=10, search="grid")###
  trmod = train(yinfect~.,data=X.use,ntree=1000,method = "rf",tuneGrid = tunegrid,trControl = trcon,metric='Accuracy')
  pred = predict(trmod,newdata = X.test,type = "prob")
  rocm <- roc(Y.test,pred[,2])
  opt.thre = as.numeric(coords(strc, "best", ret = "threshold"))
  pred.class <- ifelse(pred[,2] > opt.thre,1,0)
  dt <- table(factor(pred.class,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
  T = confusionMatrix(dt)
  return(c(T$overall[1],T$byClass[1:2]))
}

res.asc.con.rf = do.call(rbind,lapply(1:10,rate_calculate,X.filtered,final.select.asc.con))
res.asc.bin.rf = do.call(rbind,lapply(1:10,rate_calculate,X.binary.filtered,final.select.asc.bin))
res.reinf.con.rf = do.call(rbind,lapply(1:10,rate_calculate,XX.filtered,final.select.reinf.con))
res.reinf.bin.rf = do.call(rbind,lapply(1:10,rate_calculate,XX.binary.filtered,final.select.reinf.bin))
res.uninf.con.rf = do.call(rbind,lapply(1:10,rate_calculate,XXX.filtered,final.select.uninf.con))
res.uninf.bin.rf = do.call(rbind,lapply(1:10,rate_calculate,XXX.binary.filtered,final.select.uninf.bin))
