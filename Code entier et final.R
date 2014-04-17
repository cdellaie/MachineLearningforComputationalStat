library(igraph)
library(Matrix)
library(MASS)
library(kernlab)
library(ROCR)

########Importation et traitement des données

xtrain = read.table("xtrain.txt")
nomGen=xtrain[,1]
xtrain=t(as.matrix(xtrain[,-1]))
rownames(xtrain)=NULL
colnames(xtrain)=nomGen
difGen=setdiff(colnames(xtrain),colnames(lap))


ytrain = as.matrix(read.table("ytrain.txt"))

xtest = read.table("xtest.txt")
nomGen=xtest[,1]
xtest=t(as.matrix(xtest[,-1]))
rownames(xtest)=NULL
colnames(xtest)=nomGen


########### I. SVM et cross-validations

## Fonctions préliminaires

# Découpage de l'échantillon en nfolds sous échantillon de manière aléatoire
cv.folds <- function(n,nfolds=3)
{
  return(split(sample(n),rep(1:nfolds,length=n)))
}

# Prédiction de la cross-validation
cvpred.ksvm <- function(x,y,folds=3,predtype="response",...)
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

# AUC résultant de cross-validations moyenné selon n simulation
cv.auc=function(x,y,C,n=10,...)
{
  som=0
  for(i in 1:n)
  {
    print(i)
    cv=cvpred.ksvm(x,y,folds=5,predtype="response",C=C,...)
    cv.pred=prediction(cv,y)
    cv.perf=performance(cv.pred,measure="auc")
    som=cv.perf@y.values[[1]]+som
  }
  return(som/n)
}

## Cross-validations des marges C allant de 2^-10 à 2^10 selon différents noyaux

# Noyau gaussien
cv.gaus=data.frame(C=2^seq(-10,10))
cv.gaus$res=sapply(cv.gaus$C,function(C) {
  print(C)
  print(Sys.time())
  cv.auc(xtrain,ytrain,C,n=10,type="C-svc",kernel="rbfdot")
})
write.csv2(cv.gaus,"cv_gaus.csv" ,row.names=F)

# Noyau linéaire
cv.lin=data.frame(C=2^seq(-10,10))
cv.lin$res=sapply(cv.lin$C,function(C) {
  print(C)
  print(Sys.time())
  cv.auc(xtrain,ytrain,C,n=10,type="C-svc",kernel="vanilladot")
})
write.csv2(cv.lin,"cv_lin.csv" ,row.names=F)

#Noyau Laplacien
cv.lap=data.frame(C=2^seq(-10,10))
cv.lap$res=sapply(cv.lap$C,function(C) {
  print(C)
  print(Sys.time())
  cv.auc(xtrain,ytrain,C,n=10,type="C-svc",kernel="laplacedot")
})
write.csv2(cv.lap,"cv_lap.csv" ,row.names=F)

#Noyau tangente hyperbolique
cv.tanh=data.frame(C=2^seq(-10,10))
cv.tanh$res=sapply(cv.tanh$C,function(C) {
  print(C)
  print(Sys.time())
  cv.auc(xtrain,ytrain,C,n=10,type="C-svc",kernel="tanhdot")
})
write.csv2(cv.tanh,"cv_tanh.csv" ,row.names=F)

#Noyau bessélien
cv.bessel=data.frame(C=2^seq(-10,10))
cv.bessel$res=sapply(cv.bessel$C,function(C) {
  print(C)
  print(Sys.time())
  cv.auc(xtrain,ytrain,C,n=10,type="C-svc",kernel="besseldot")
})
write.csv2(cv.bessel,"cv_bessel.csv" ,row.names=F)

# Noyau polynomiale de degré 2
cv.poly=data.frame(C=2^seq(-10,10))
cv.poly$res=sapply(cv.poly$C,function(C) {
  print(C)
  print(Sys.time())
  cv.auc(xtrain,ytrain,C,type="C-svc",kernel=polydot(degree=2),scaled=c())
})
write.csv2(cv.poly,"cv.poly.csv" ,row.names=F)

# Noyau Anova
cv.anova=data.frame(C=2^seq(-10,10))
cv.anova$res=sapply(cv.anova$C,function(C) {
  print(C)
  print(Sys.time())
  cv.auc(xtrain,ytrain,C,n=1,type="C-svc",kernel="anovadot",sigma=1,degre=1,scaled=c())
})
write.csv2(cv.anova,"cv_anova.csv" ,row.names=F)


########### II. Utilisation des intéractions entre gènes comme a-priori

##Importation et traitement des données

ppi=read.table("ppi.txt")
ppi.reduce=ppi[(ppi[,1] %in% nomGen)&(ppi[,2] %in% nomGen),]
genes=unique(c(ppi.reduce[,1],ppi.reduce[,2]))

nomGenRe=c(colnames(lap),setdiff(nomGen,colnames(lap)))
xtrainRe=xtrain[,nomGenRe]

### Calcul du Laplacien et de ses inverses

## Laplacien
greduc=graph.data.frame(ppi.reduce)
lap=graph.laplacian(greduc,sparse=F)
n=dim(lap)[1]

## Inverse avec B=0

tes=solve(lapS[1:n,1:n]+diag(1/n^2,n))
tes2=diag(0,dim(xtrain)[2])
tes2[1:n,1:n]=as.matrix(tes)

# Noyau associé
kerMat=xtrainRe%*%tes2%*%t(xtrainRe)
write.table(kerMat,"kerMat.txt")

## Inverses et noyaux associés avec B allant de 10^-3 à 10^3

for (C in 10^seq(-3,3))
{
  print(C)
  lapsC=diag(C,dim(xtrain)[2])
  lapsC[1:n,1:n]=lap+diag(C,n)
  tes2=solve(lapsC)
  kerMat=xtrainRe%*%tes2%*%t(xtrainRe)
  write.table(kerMat,paste("kerMat",C,".txt",sep=""))
  
}

### Cross-validations des noyaux

## Fonctions préléminaires

cvpred.ksvmK <- function(x,y,folds=3,predtype="response",...)
{
  n <- length(y)
  ypred <- numeric(n)
  s <- cv.folds(n,folds)
  for (i in seq(folds)) {
    m <- ksvm(as.kernelMatrix(x[-s[[i]],-s[[i]]]),y[-s[[i]]],...)
    ypred[s[[i]]] <- predict(m,as.kernelMatrix(x[s[[i]],-s[[i]]]),type=predtype)
  }
  invisible(ypred)
}

cv.aucK=function(x,y,C,n=10,...)
{
  som=0
  for(i in 1:n)
  {
    print(i)
    cv=cvpred.ksvmK(x,y,folds=5,predtype="response",C=C,...)
    cv.pred=prediction(cv,y)
    cv.perf=performance(cv.pred,measure="auc")
    som=cv.perf@y.values[[1]]+som
  }
  return(som/n)
}

## Cross-validations des marges C allant de 2^-10 à 2^10 selon différents noyaux

# B=0
kermat=as.matrix(read.table("kermat.txt"))
dimnames(kermat)=NULL
cv.ker0=data.frame(C=2^seq(-10,10))
cv.ker0$res=sapply(cv.ker0$C,function(C) {
  print(C)
  print(Sys.time())
  cv.auc(kermat,ytrain,C,n=1,type="C-svc")
})
write.csv2(cv.ker0,"cv_ker0.csv" ,row.names=F)

# B allant de 10^-3 à 10^3
for (i in 10^seq(-3,3))
{
  print(paste("kermat",i,".txt"))
  kermat=as.matrix(read.table(paste("kermat",i,".txt",sep="")))
  dimnames(kermat)=NULL
  cv.ker0=data.frame(C=2^seq(-10,10))
  cv.ker0$res=sapply(cv.ker0$C,function(C) {
    print(C)
    print(Sys.time())
    cv.auc(kermat,ytrain,C,n=1,type="C-svc")
  })
  write.csv2(cv.ker0,paste("cv_ker",i,".csv",sep="") ,row.names=F)
}

########### III. Conclusion avec les prédictions finales

## Modèle sélectionné : noyau polynomiale de degré 2 avec une marge C=64
modFin=ksvm(xtrain,ytrain,C=64,type="C-svc",kernel=polydot(degree=2),scaled=c())

##Prédictions sur les 92 individus tests
finRes=predict(modFin,xtest,type="decision")
write.table(finRes,"prédictions.txt",col.names=F,row.names=F)
