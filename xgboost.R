library(caret)
library(xgboost)

xgb.varimp <- function(j,X1)   ####calculate variable importance
{
  set.seed(232+7*j)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.train = data.frame(X1[train_ind,-1])
  Y.train = as.factor(X1[train_ind,1])
  xgtrain <- xgb.DMatrix(data = as.matrix(X.train),label = as.numeric(Y.train)-1)
  xgb1 = xgb.train(data = xgtrain, booster = "gblinear",nrounds=500, watchlist=list(data = xgtrain), objective = "binary:logistic",verbose = FALSE)
  xgi <- xgb.importance(model = xgb1)
  #trcon = trainControl(method="cv", number=10, search="grid")
  #xgmod = train(X.train,Y.train,method="xgbLinear",trControl = trcon)
  #xgs = varImp(xgmod)$importance
  #x = xgs[,1]
  xgi = xgi[match(colnames(X.train),xgi$Feature),]
  x = xgi$Weight
  names(x) = xgi$Feature
  return(x)
}

params <- list(booster = "gblinear", 
               objective = "binary:logistic")
xgb_base <- xgb.train (params = params,
                       data = xgd,
                       nrounds =1000,
                       print_every_n = 10,
                       eval_metric = "auc",
                       eval_metric = "error",
                       early_stopping_rounds = 50)

xg.var.rank.asc.con = sort(abs(rowMeans(do.call(cbind,lapply(1:10,xgb.varimp,X.filtered)))),decreasing = TRUE)
xg.var.rank.asc.bin = sort(abs(rowMeans(do.call(cbind,lapply(1:10,xgb.varimp,X.binary.filtered)))),decreasing = TRUE)
xg.var.rank.reinf.con = sort(abs(rowMeans(do.call(cbind,lapply(1:10,xgb.varimp,XX.filtered)))),decreasing = TRUE)
xg.var.rank.reinf.bin = sort(rowMeans(do.call(cbind,lapply(1:10,xgb.varimp,XX.binary.filtered))),decreasing = TRUE)
xg.var.rank.uninf.con = sort(abs(rowMeans(do.call(cbind,lapply(1:10,xgb.varimp,XXX.filtered)))),decreasing = TRUE)
xg.var.rank.uninf.bin = sort(rowMeans(do.call(cbind,lapply(1:10,xgb.varimp,XXX.binary.filtered))),decreasing = TRUE)

num.calculate.xg <- function(i,var.rank,X1) ### calculate best number in the final model
{
  select.var <- names(var.rank)[1:i]
  error = 0
  for(j in 1:10)
  {
    set.seed(232+7*j)
    train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
    X.train.use <- X1[train_ind,select.var]
    Y.train.use = as.factor(X1[train_ind,1])
    xgtrain <- xgb.DMatrix(data = as.matrix(X.train.use),label = as.numeric(Y.train.use)-1)
    #trcon = trainControl(method="repeatedcv", number=5,repeats = 3 ,search="grid")
    xgb1 = xgb.cv(data = xgtrain,eval_metric = "error",
                  early_stopping_rounds= 50,nfold = 10, booster = "gblinear",nrounds=500, watchlist=list(data = xgtrain), objective = "binary:logistic",verbose = FALSE)

    #nbmod = train(X.train.use,Y.train.use,method="xgbLinear",trControl = trcon,tuneGrid = grid)
    error[j] <- min(xgb1$evaluation_log$test_error_mean)
    #max(nbmod$results[,"Accuracy"],na.rm=TRUE)
  }
  perform <- mean(error,na.rm=TRUE)
  return(perform)
  #tunegrid <- expand.grid(.mtry=1:i)
  #trcon = trainControl(method="repeatedcv", number=10, repeats=3,search="grid")
  #trmod = train(yinfect~.,data=X.use,ntree=1000,tuneGrid = tunegrid,trControl = trcon,metric='Accuracy')
}

# grid <- expand.grid(
#   nrounds =1000,
#   eta = c(0.1, 0.3),
#   lambda = c(0.5,1),
#   alpha = c(0,0.5)
# )
xg.perform.res.asc.con = unlist(lapply(2:35,num.calculate.xg,xg.var.rank.asc.con,X.filtered))

xg.perform.res.asc.bin = unlist(lapply(2:35,num.calculate.xg,xg.var.rank.asc.bin,X.binary.filtered))
xg.perform.res.reinf.con = unlist(lapply(2:35,num.calculate.xg,xg.var.rank.reinf.con,XX.filtered))

xg.perform.res.reinf.bin = unlist(lapply(2:35,num.calculate.xg,xg.var.rank.reinf.bin,XX.binary.filtered))

xg.perform.res.uninf.con = unlist(lapply(2:35,num.calculate.xg,xg.var.rank.uninf.con,XXX.filtered))

xg.perform.res.uninf.bin = unlist(lapply(2:35,num.calculate.xg,xg.var.rank.uninf.bin,XXX.binary.filtered))

xg.final.select.asc.con <- names(xg.var.rank.asc.con)[1:(which.min(xg.perform.res.asc.con)+1)]
xg.final.select.asc.bin <- names(xg.var.rank.asc.bin)[1:(which.min(xg.perform.res.asc.con)+1)]
                                                         
rate_calculate.xgb<-  function(j,select.var,X1)
  {
    
    set.seed(232+7*j)
    train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
    X.test <- X1[-train_ind,select.var]
    Y.test = as.factor(X1[-train_ind,1])
    X.train = data.frame(X1[train_ind,-1])
    Y.train = as.factor(X1[train_ind,1])
    X.use <- X1[train_ind,c("yinfect",select.var)]
    #grid <- expand.grid(fL=c(0,0.5,1.0), usekernel = TRUE, adjust=c(0,0.5,1.0))
    grid <- expand.grid(
      nrounds =500,
      eta =  0.3,
      lambda = 1,
      alpha = c(0,0.5,1)
    )
    #xgtrain <- xgb.DMatrix(data = as.matrix(X.train),label = as.numeric(Y.train)-1)
    #xgtest <- xgb.DMatrix(data = as.matrix(X.test),label = as.numeric(Y.test)-1)
    # xgb1 = xgb.train(data = xgtrain, booster = "gblinear",nrounds=500, watchlist=list(train = xgtrain,eval = xgtest), objective = "binary:logistic",verbose = FALSE)
    #pred = predict(xgb1,newdata = xgtest,type="prob")
    X.use <- X.use %>% mutate(yinfect = factor(yinfect, 
                                               labels = make.names(levels(yinfect))))
    trcon = trainControl(classProbs = TRUE,method="cv", number=10,search="grid",summaryFunction = twoClassSummary)
    trmod = train(yinfect~.,data=X.use,method="xgbLinear",trControl = trcon,metric='ROC',tuneGrid = grid)
    pred = predict(trmod,newdata = X.test,type="prob")
    strc <- roc(Y.test,pred[,2])
    opt.thre = as.numeric(coords(strc, "best", ret = "threshold",transpose = TRUE))
    pred.class <- ifelse(pred[,2] > opt.thre,1,0)
    dt <- table(factor(pred.class,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
    T = confusionMatrix(dt)
    return(c(T$overall[1],T$byClass[1:2]))
}

res.asc.con.xg = do.call(rbind,lapply(1:10,rate_calculate.xgb,xg.final.select.asc.con,X.filtered))
res.asc.bin.xg = do.call(rbind,lapply(1:10,rate_calculate.xgb,xg.final.select.asc.bin,X.binary.filtered))
res.reinf.con.xg = do.call(rbind,lapply(1:10,rate_calculate.xgb,xg.final.select.reinf.con,XX.filtered))
res.reinf.bin.xg = do.call(rbind,lapply(1:10,rate_calculate.xgb,xg.final.select.reinf.bin,XX.binary.filtered))
res.uninf.con.xg = do.call(rbind,lapply(1:10,rate_calculate.xgb,xg.final.select.uninf.con,XXX.filtered))                                                        
res.uninf.bin.xg = do.call(rbind,lapply(1:10,rate_calculate.xgb,xg.final.select.uninf.bin,XXX.binary.filtered))