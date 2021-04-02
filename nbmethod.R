library(caret)
library(tidyverse)

X1 = na.omit(X.binary)
X1[,1]=as.factor(X1[,1])
XXX <- na.omit(bin.uninfect)
id <- apply(XXX,2,function(x){sum(x!=0,na.rm=TRUE)})

XXX.binary.filtered <- XXX[,names(id[id>2])]


nb.varimp <- function(j,X1)   ####calculate variable importance
{
  set.seed(232+7*j)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.train = data.frame(X1[train_ind,-1])
  Y.train = as.factor(X1[train_ind,1])
  trcon = trainControl(method="cv", number=5,search="grid")
  grid <- expand.grid(fL=c(0,0.5,1.0), usekernel = TRUE, adjust=c(0,0.5,1.0))
  nbmod = train(X.train,Y.train,method="nb",trControl = trcon)#tuneGrid = grid)
  nbs = varImp(nbmod)$importance
  x = nbs[,1]
  names(x) = rownames(nbs)
  return(x)
}

nb.var.rank.asc.con = sort(rowMeans(do.call(cbind,lapply(1:10,nb.varimp,X.filtered))),decreasing = TRUE)
nb.var.rank.asc.bin= sort(rowMeans(do.call(cbind,lapply(1:10,nb.varimp,X.binary.filtered))),decreasing = TRUE)

nb.var.rank.reinf.con = sort(rowMeans(do.call(cbind,lapply(1:10,nb.varimp,XX.filtered))),decreasing = TRUE)
nb.var.rank.reinf.bin = sort(rowMeans(do.call(cbind,lapply(1:10,nb.varimp,XX.binary.filtered))),decreasing = TRUE)

nb.var.rank.uninf.bin = sort(rowMeans(do.call(cbind,lapply(1:10,nb.varimp,XXX.binary.filtered))),decreasing = TRUE)
nb.var.rank.uninf.con = sort(rowMeans(do.call(cbind,lapply(1:10,nb.varimp,XXX.filtered))),decreasing = TRUE)

num.calculate.nb <- function(i,var.rank,X1) ### calculate best number in the final model
{
  select.var <- names(var.rank)[1:i]
  accu = 0
  for(j in 1:10)
  {
    set.seed(232+7*j)
    train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
    X.train.use <- X1[train_ind,select.var]
    Y.train.use = as.factor(X1[train_ind,1])
    X.use <- X1[train_ind,c("yinfect",select.var)]
    trcon = trainControl(method="cv", number=10,search="grid",classProbs = TRUE,summaryFunction = twoClassSummary)
    grid <- expand.grid(fL=0, usekernel = TRUE, adjust=c(0.5,1.0))
    X.use <- X.use %>% mutate(yinfect = factor(yinfect, 
                                                  labels = make.names(levels(yinfect))))
    nbmod = train(yinfect~.,data=X.use,method="nb",trControl = trcon,tuneGrid = grid,metric='ROC')
    accu[j] <- max(nbmod$results[,"ROC"],na.rm=TRUE)
  }
  perfom <- mean(accu,na.rm=TRUE)
  #tunegrid <- expand.grid(.mtry=1:i)
  #trcon = trainControl(method="repeatedcv", number=10, repeats=3,search="grid")
  #trmod = train(yinfect~.,data=X.use,ntree=1000,tuneGrid = tunegrid,trControl = trcon,metric='Accuracy')
}

nb.perform.res.asc.con = unlist(lapply(2:35,num.calculate.nb,nb.var.rank.asc.con,X.filtered))
nb.perform.res.asc.bin = unlist(lapply(2:35,num.calculate.nb,nb.var.rank.asc.bin,X.binary.filtered))

nb.perform.res.reinf.con = unlist(lapply(2:35,num.calculate.nb,nb.var.rank.reinf.con,XX.filtered))
nb.perform.res.reinf.bin = unlist(lapply(2:35,num.calculate.nb,nb.var.rank.reinf.bin,XX.binary.filtered))
nb.perform.res.uninf.con = unlist(lapply(2:35,num.calculate.nb,nb.var.rank.uninf.con,XXX.filtered))

nb.perform.res.uninf.bin = unlist(lapply(2:35,num.calculate.nb,nb.var.rank.uninf.bin,XXX.binary.filtered))

nb.final.select.asc.con <- names(nb.var.rank.asc.con)[1:(which.max(nb.perform.res.asc.con)+1)]
nb.final.select.asc.bin <- names(nb.var.rank.asc.bin)[1:(which.max(nb.perform.res.asc.bin)+1)]

nb.final.select.reinf.con <- names(nb.var.rank.reinf.con)[1:(which.max(nb.perform.res.reinf.con)+1)]
nb.final.select.reinf.bin <- names(nb.var.rank.reinf.bin)[1:(which.max(nb.perform.res.reinf.bin)+1)]
nb.final.select.uninf.bin <- names(nb.var.rank.uninf.bin)[1:(which.max(nb.perform.res.uninf.bin)+1)]
nb.final.select.uninf.con <- names(nb.var.rank.uninf.con)[1:(which.max(nb.perform.res.uninf.con)+1)]

rate_calculate.nb <- function(j,select.var,X1)
{
  
  set.seed(232+7*j)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.test <- X1[-train_ind,select.var]
  Y.test = as.factor(X1[-train_ind,1])
  X.use <- X1[train_ind,c("yinfect",select.var)]
  X.use <- X.use %>% mutate(yinfect = factor(yinfect, 
                                             labels = make.names(levels(yinfect))))
  
  grid <- expand.grid(fL=0, usekernel = TRUE, adjust=c(0.5,1.0))
  trcon = trainControl(classProbs = TRUE,method="cv",number = 10, search="grid",summaryFunction = twoClassSummary)
  trmod = train(yinfect~.,data=X.use,method="nb",trControl = trcon,metric='ROC',tuneGrid = grid)
  pred = predict(trmod,newdata = X.test,type="prob")
  rocm <- roc(Y.test,pred[,2])
  opt.thre = as.numeric(coords(rocm, "best", ret = "threshold",transpose = TRUE))
  pred.class <- ifelse(pred[,2] > opt.thre,1,0)
  dt <- table(factor(pred.class,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
  T = confusionMatrix(dt)
  return(c(T$overall[1],T$byClass[1:2]))
}

res.asc.con.nb = do.call(rbind,lapply(1:10,rate_calculate.nb,nb.final.select.asc.con,X.filtered))
res.asc.bin.nb = do.call(rbind,lapply(1:10,rate_calculate.nb,nb.final.select.asc.bin,X.binary.filtered))

res.reinf.con.nb = do.call(rbind,lapply(1:10,rate_calculate.nb,nb.final.select.reinf.con,XX.filtered))
res.reinf.bin.nb = do.call(rbind,lapply(1:10,rate_calculate.nb,nb.final.select.reinf.bin,XX.binary.filtered))

res.uninf.con.nb1= do.call(rbind,lapply(1:10,rate_calculate.nb,nb.final.select.uninf.con,XXX.filtered))

res.uninf.bin.nb = do.call(rbind,lapply(1:10,rate_calculate.nb,nb.final.select.uninf.bin,XXX.binary.filtered))

continuous_res_asc <- data.frame(rbind(res.asc.con.rf,res.asc.con.nb,res.asc.con.xg,res.asc.con.svm,res.asc.con.knn,res.asc.con.els),
                                 method=c(rep("Random Forest",10),rep("Naive Bayes",10),rep("xgbLinear",10),rep("msvm-rfe",10),rep("knn",10),rep("Elastic Net",10)))
dat.asc.con <- reshape2::melt(continuous_res_asc)

continuous_res_reinf <- data.frame(rbind(res.reinf.con.rf,res.reinf.con.nb,res.reinf.con.xg,res.reinf.con.svm,res.reinf.con.knn,res.reinf.con.els),
                                   method=c(rep("Random Forest",10),rep("Naive Bayes",10),rep("xgbLinear",10),rep("msvm-rfe",10),rep("knn",10),rep("Elastic Net",10)))
dat.reinf.con <- reshape2::melt(continuous_res_reinf)

continuous_res_uninf <- data.frame(rbind(res.uninf.con.rf,res.uninf.con.nb,res.uninf.con.xg,res.uninf.con.svm,res.uninf.con.knn,res.uninf.con.els1),
                                   method=c(rep("Random Forest",10),rep("Naive Bayes",10),rep("xgbLinear",10),rep("msvm-rfe",10),rep("knn",10),rep("Elastic Net",10)))
dat.uninf.con <- reshape2::melt(continuous_res_uninf)

binary_res_uninf <- data.frame(rbind(res.uninf.bin.rf,res.uninf.bin.nb,res.uninf.bin.xg,res.uninf.bin.svm,res.uninf.bin.knn,res.uninf.bin.els),
                               method=c(rep("Random Forest",10),rep("Naive Bayes",10),rep("xgbLinear",10),rep("msvm-rfe",10),rep("knn",10),rep("Elastic Net",10)))
dat.uninf.bin <- reshape2::melt(binary_res_uninf)
binary_res_reinf <- data.frame(rbind(res.reinf.bin.rf,res.reinf.bin.nb,res.reinf.bin.xg,res.reinf.bin.svm,res.reinf.bin.knn,res.reinf.bin.els),
                               method=c(rep("Random Forest",10),rep("Naive Bayes",10),rep("xgbLinear",10),rep("msvm-rfe",10),rep("knn",10),rep("Elastic Net",10)))
dat.reinf.bin <- reshape2::melt(binary_res_reinf)

binary_res_asc <- data.frame(rbind(res.asc.bin.rf,res.asc.bin.nb,res.asc.bin.xg,res.asc.bin.svm,res.asc.bin.knn,res.asc.bin.els),method=c(rep("Random Forest",10),rep("Naive Bayes",10),rep("xgbLinear",10),rep("msvm-rfe",10),rep("knn",10),rep("Elastic Net",10)))
dat.asc.bin <- reshape2::melt(binary_res_asc)
ggplot(dat.asc.bin,aes(x = method, y =value, fill = method))+geom_boxplot()+scale_x_discrete(breaks = NULL)+
  facet_grid(cols = vars(variable))+labs(title= "ascender binary",x="",y="")+theme_bw()+
  theme(axis.text.y=element_text(size=15,face = "bold"),axis.ticks = element_blank(), axis.text.x = element_blank(),strip.text = element_text(size = 20))


ggplot(dat.asc.con,aes(x = method, y =value, fill = method))+geom_boxplot()+scale_x_discrete(breaks = NULL)+
  facet_grid(cols = vars(variable))+labs(title= "ascender continuous",x="",y="")+theme_bw()+
  theme(axis.text.y=element_text(size=15,face = "bold"),axis.ticks = element_blank(), axis.text.x = element_blank(),strip.text = element_text(size = 20))

ggplot(dat.asc.bin,aes(x = method, y =value, fill = method))+geom_boxplot()+scale_x_discrete(breaks = NULL)+
  facet_grid(cols = vars(variable))+labs(title= "ascender binary",x="",y="")+theme_bw()+
  theme(axis.text.y=element_text(size=15,face = "bold"),axis.ticks = element_blank(), axis.text.x = element_blank(),strip.text = element_text(size = 20))


ggplot(dat.reinf.con,aes(x = method, y =value, fill = method))+geom_boxplot()+scale_x_discrete(breaks = NULL)+
  facet_grid(cols = vars(variable))+labs(title= "reinfect continuous",x="",y="")+theme_bw()+
  theme(axis.text.y=element_text(size=15,face = "bold"),axis.ticks = element_blank(), axis.text.x = element_blank(),strip.text = element_text(size = 20))
ggplot(dat.reinf.bin,aes(x = method, y =value, fill = method))+geom_boxplot()+scale_x_discrete(breaks = NULL)+
  facet_grid(cols = vars(variable))+labs(title= "reinfect binary",x="",y="")+theme_bw()+
  theme(axis.text.y=element_text(size=15,face = "bold"),axis.ticks = element_blank(), axis.text.x = element_blank(),strip.text = element_text(size = 20))


ggplot(dat.asc.con,aes(x = method, y =value, fill = method))+geom_boxplot()+scale_x_discrete(breaks = NULL)+
  facet_grid(cols = vars(variable))+labs(title= "ascender continuous",x="",y="")+theme_bw()+
  theme(axis.text.y=element_text(size=15,face = "bold"),axis.ticks = element_blank(), axis.text.x = element_blank(),strip.text = element_text(size = 20))

ggplot(dat.uninf.con,aes(x = method, y =value, fill = method))+geom_boxplot()+scale_x_discrete(breaks = NULL)+
  facet_grid(cols = vars(variable))+labs(title= "reinfect-Uninfected continuous",x="",y="")+theme_bw()+
  theme(axis.text.y=element_text(size=15,face = "bold"),axis.ticks = element_blank(), axis.text.x = element_blank(),strip.text = element_text(size = 20))
ggplot(dat.uninf.bin,aes(x = method, y =value, fill = method))+geom_boxplot()+scale_x_discrete(breaks = NULL)+
  facet_grid(cols = vars(variable))+labs(title= "reinfect-Uninfected binary",x="",y="")+theme_bw()+
  theme(axis.text.y=element_text(size=15,face = "bold"),axis.ticks = element_blank(), axis.text.x = element_blank(),strip.text = element_text(size = 20))

binary_res_reinf <- data.frame(rbind(res.reinf.bin.rf,res.reinf.bin.nb,res.reinf.bin.xg,res.reinf.bin.svm),method=c(rep("rf",10),rep("nb",10),rep("xgbLinear",10)))
dat.reinf.bin <- reshape2::melt(binary_res_reinf)
ggplot(dat.reinf.bin,aes(x = method, y =value, fill = method))+geom_boxplot()+
  facet_grid(cols = vars(variable))+labs(title= "reinfect binary")



#####Top 5 boxplot######

topbox <- function(select,X)
{
  se <- select[1:min(5,length(select))]
  dat.use <- X[,c("yinfect",se)]
  dtu <- reshape2::melt(dat.use,id.vars = "yinfect")
  if(grepl(pattern = "asc",substitute(select)))
  {
    p = ggplot(dtu,aes(x = yinfect,y = value,col=yinfect))+
      geom_boxplot()+geom_jitter(width = 0.4,height = 0.05)+
      facet_wrap(~variable,ncol = 3,scales = "free")+
      scale_color_manual(values = c("#00BFC4","#F8766D"),name = "Infection Status",labels = c("Endo-","Endo+"))+labs(x="",y="")+theme(strip.text = element_text(size = 18))
    
  }
  if(grepl(pattern = "reinf",substitute(select)))
  {
    p = ggplot(dtu,aes(x = yinfect,y = value,col=yinfect))+geom_boxplot()+geom_jitter(width = 0.4,height = 0.05)+
      facet_wrap(~variable,ncol = 3,scales = "free")+
      scale_color_manual(values = c("#00BFC4","#F8766D"),name = "Infection Status",labels = c("Protective","Susceptible"))+labs(x="",y="")+theme(strip.text = element_text(size = 18))
  }
  else if(grepl(pattern = "uninf",substitute(select)))
  {
    p = ggplot(dtu,aes(x = yinfect,y = value,col=yinfect))+geom_boxplot()+geom_jitter(width = 0.4,height = 0.05)+
      facet_wrap(~variable,ncol = 3,scales = "free")+
      scale_color_manual(values = c("#00BFC4","#F8766D"),name = "Infection Status",labels = c("Protective_Uninf","Susceptible_Uninf"))+labs(x="",y="")+theme(strip.text = element_text(size = 18))
  }
  return(p)
}


cifunc <- function(x,i)
{
  m = colMeans(x)
  st <- apply(x,2,sd)
  lower = m[i] - 1.96*st[i]
  upper = m[i] + 1.96*st[i]
  return(c(m[i],lower,upper))
}