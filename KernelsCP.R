.libPaths("C:/Users/Clément/Documents/R-library")
setwd("C:/Users/Clément/Documents/GitHub/MachineLearningforComputationalStat")
# install.packages("kernlab",lib="C:/Users/Clément/Documents/R-library")
# install.packages("ROCR",lib="C:/Users/Clément/Documents/R-library")
# install.packages("igraph",lib="C:/Users/Clément/Documents/R-library")
library(kernlab)
library(ROCR)
library(igraph)

#téléchargement des données

xtrain = read.table("xtrain.txt")
nomGen=xtrain[,1]
xtrain=t(as.matrix(xtrain[,-1]))
rownames(xtrain)=NULL
colnames(xtrain)=nomGen

ytrain = as.matrix(read.table("ytrain.txt"))

test = read.table("xtest.txt")
nomGen=test[,1]
xtest=t(as.matrix(test[,-1]))
rownames(xtest)=NULL
colnames(xtest)=nomGen

cv.folds <- function(n,nfolds=3)
  ## Randomly split the n samples into folds
  ## Returns a list of nfolds lists of indices, each corresponding to a fold
{
  return(split(sample(n),rep(1:nfolds,length=n)))
}

cvpred.ksvm <- function(x,y,folds=3,predtype="response",...)
  ## Return a vector of predictions by cross-validation
  ## 'predtype' should be one of response (by default), decision or probabilities, depending the prediction we want (SVM label, score or probability, see predict.ksvm())
  ## Additional parameters are passed to ksvm() to train the SVM
{
  n <- length(y)
  ypred <- numeric(n)
  s <- cv.folds(n,folds)
  for (i in seq(folds)) {
    m <- ksvm(x[-s[[i]],],y[-s[[i]]],...)
    ypred[s[[i]]] <- predict(m,x[s[[i]],],type=predtype)
  }
  invisible(ypred)
}

cv.auc=function(x,y,Cte,...)
{
  cv=cvpred.ksvm(x,y,folds=5,predtype="response",C=Cte,...)
  cv.pred=prediction(cv,y)
  cv.perf=performance(cv.pred,measure="auc")
  return(cv.perf@y.values[[1]])
}

###Cross validation pour C allant de 2^-10 à 2^10
cv.poly=data.frame(C=2^seq(-10,10))

cv.auc(datat,ytrain,C=1,type="C-svc",kernel=polydot(degree=2),scaled=c())
ksvm(datat,ytrain,type="C-svc",kernel=polydot(degree=2),scale=1)


cv.poly$res=sapply(cv.poly$C,function(C) cv.auc(xtrain,ytrain,C,type="C-svc",kernel=polydot(degree=2),scaled=c()))
par(xlog=T)
plot(cv.poly$C,cv.poly$res,type="l",log="x")

##Meilleur C trouvé
Cmax=cv.poly$C[which.max(cv.poly$res)]



#à virer  xtrain[,1]




