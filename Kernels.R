library(kernlab)
library(ROCR)
library(igraph)

##########################################
#téléchargement et traitement des données#
##########################################

xtrain = read.table("C:/Users/Clément/Desktop/ENSAE 3A/MLCB/xtrain.txt")
ytrain = read.table("C:/Users/Clément/Desktop/ENSAE 3A/MLCB/ytrain.txt")
test = read.table("C:/Users/Clément/Desktop/ENSAE 3A/MLCB/xtest.txt")

matrix(test)
ttest=as.matrix(test[2:93])
FinalTest=t(ttest)

matrix(xtrain)
data=as.matrix(xtrain[,2:185])
datat=t(data)

#creation d un set d'entraînement et d'un set test
xtrain2=datat[1:150,]
ytrain2=ytrain[1:150,]

xtest=datat[151:184,]
ytest=ytrain[151:184,]

##################################################################################################################
#essai sur les set d'entraînements avec plusieurs noyaux                                                        ##
#le meilleur sera choisi pour être entrainé sur la base complète et prévoir sur les données à traiter           ## 
#tanhdot,laplacedot,besseldot,anovadot                                                                          ##
##################################################################################################################

LinearKernel <- ksvm(xtrain2,ytrain2,type="C-svc",kernel="vanilladot")
plot(LinearKernel,data=xtest)
ypred=predict(LinearKernel,xtest)
table(ytest,ypred)
sum(ypred==ytest)/length(ytest)

GaussianKernel <- ksvm(xtrain2,ytrain2,type="C-svc",kernel="rbfdot",kpar=list(sigma=1),scale=FALSE)
plot(GaussianKernel,data=xtest)
ypred=predict(GaussianKernel,xtest)
table(ytest,ypred)
sum(ypred==ytest)/length(ytest)

PolyKernel = ksvm(xtrain2,ytrain2,type="C-svc",kernel=polydot(degree=2))
plot(PolyKernel,data=xtest)
ypred=predict(PolyKernel,xtest)
table(ytest,ypred)
sum(ypred==ytest)/length(ytest)

AnovaKernel = ksvm(xtrain2,ytrain2,type="C-svc",kernel="anovadot")
plot(AnovaKernel,data=xtest)
ypred=predict(AnovaKernel,xtest)
table(ytest,ypred)
sum(ypred==ytest)/length(ytest)

TanhKernel = ksvm(xtrain2,ytrain2,type="C-svc",kernel="tanhdot")
plot(TanhKernel,data=xtest)
ypred=predict(TanhKernel,xtest)
table(ytest,ypred)
sum(ypred==ytest)/length(ytest)

BesselKernel = ksvm(xtrain2,ytrain2,type="C-svc",kernel="besseldot")
plot(BesselKernel,data=xtest)
ypred=predict(BesselKernel,xtest)
table(ytest,ypred)
sum(ypred==ytest)/length(ytest)

LaplaceKernel = ksvm(xtrain2,ytrain2,type="C-svc",kernel="laplacedot")
plot(LaplaceKernel,data=xtest)
ypred=predict(LaplaceKernel,xtest)
table(ytest,ypred)
sum(ypred==ytest)/length(ytest)

#############################################################################
#Linear, Anova et Poly donnent les meilleurs résultats : 0.67 accuracy ######
#############################################################################

BestKernel = ksvm(datat,ytrain,type="C-svc",kernel=polydot(degree=2),scale=1)
predict(BestKernel,FinalTest)

##################################################################################################################
##################################################################################################################
###                                         Cross Validation                                                   ###
##################################################################################################################
##################################################################################################################

#Plutôt que de tester les paramètres sur un échantillon produit à la main, on se sert de la cross-validation

#############################
## Fonctions préliminaires ##
#############################

################################################################################################# 
##   Le code des fonctions préliminaires est disponible sur la page personnelle de J-P Vert    ##
##   http://cbio.ensmp.fr/~jvert                                                               ##                  
################################################################################################# 
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

#auc retourne un score de performance de l'algortihme

cv.auc=function(x,y,Cte,...)
{
  cv=cvpred.ksvm(x,y,folds=5,predtype="response",C=Cte,...)
  cv.pred=prediction(cv,y)
  cv.perf=performance(cv.pred,measure="auc")
  return(cv.perf@y.values[[1]])
}

#######################
#######################
#######################

PolyKernel = ksvm(xtrain2,ytrain2,type="C-svc",kernel=polydot(degree=2))
ypred=predict(PolyKernel,xtest)
table(ytest,ypred)
sum(ypred==ytest)/length(ytest)

ypredscore = predict(PolyKernel,xtest,type="decision")
table(ypredscore > 0,ypred)

pred <- prediction(ypredscore,ytest)
perf <- performance(pred, measure = "prec", x.measure = "rec") 
plot(perf)

# Exemple de prediction avec cross-validation

k=5	
ypredscorecv <- cvpred.ksvm(x,y,folds=k,type="C-svc",kernel=polydot(degree=2),C=1,scaled=c(),predtype="decision")

# Check the performance
print(table(ypredscorecv > 0,y))
pred <- prediction(ypredscorecv,y)

#Courbe ROC
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf)

# Estimate the CV error with ksvm directly, and compare
svp <- ksvm(x,y,type="C-svc",kernel='vanilladot',C=1,scaled=c(),cross=5)
print(cross(svp))
print(1-sum((ypredscorecv>0)==(y==1))/n)

#############################################################################
###Cross validation pour C allant de 2^-10 à 2^10                          ##
#############################################################################

xtrain3=datat[1:184,]
ytrain3=ytrain[1:184,]

###Polynôme de degré 2
cv.poly=data.frame(C=2^seq(-10,10))
cv.auc(xtrain3,ytrain3,C=1,type="C-svc",kernel=polydot(degree=2),scaled=c())
cv.poly$res=sapply(cv.poly$C,function(C) cv.auc(xtrain3,ytrain3,C,type="C-svc",kernel=polydot(degree=2),scaled=c()))
par(xlog=T)
plot(cv.poly$C,cv.poly$res,type="l",log="x")

##Meilleur C trouvé = 0.125
Cmax=cv.poly$C[which.max(cv.poly$res)]
BestPolyKernel = ksvm(datat,ytrain,type="C-svc",C=Cmax,kernel=polydot(degree=2),scale=c())
PolyPred=predict(BestPolyKernel,FinalTest)

##########
##########

###Anova
cv.auc(xtrain3,ytrain3,C=1,type="C-svc",kernel="anovadot",sigma=1,degre=1,scaled=c())
cv.poly$res=sapply(cv.poly$C,function(C) cv.auc(xtrain3,ytrain3,C,type="C-svc",kernel="anovadot",scaled=c()))
par(xlog=T)
plot(cv.poly$C,cv.poly$res,type="l",log="x")

##Meilleur C trouvé
Cmax=cv.poly$C[which.max(cv.poly$res)]
BestAnovaKernel = ksvm(datat,ytrain,type="C-svc",C=Cmax,kernel="anovadot",scale=c())
predict(BestAnovaKernel,FinalTest)

##########
##########

###Noyau gaussien
cv.auc(xtrain3,ytrain3,C=1,type="C-svc",kernel="rbfdot",kpar=list(sigma=1),scale=FALSE)
cv.poly$res=sapply(cv.poly$C,function(C) cv.auc(xtrain3,ytrain3,C,type="C-svc",kernel="rbfdot",kpar=list(sigma=1),scale=FALSE))
par(xlog=T)
plot(cv.poly$C,cv.poly$res,type="l",log="x")

##Meilleur C trouvé = 0.0009765625
Cmax=cv.poly$C[which.max(cv.poly$res)]
BestGaussianKernel = ksvm(datat,ytrain,type="C-svc",C=Cmax,kernel="rbfdot",kpar=list(sigma=1),scale=FALSE)

GaussianPred=predict(BestGaussianKernel,FinalTest)

## Estimation du meilleur margin C par moyennage (N simulations)
###Polynôme de degré 2
N=100
CC=0
cv.poly=data.frame(C=2^seq(-10,10))

for (i in 1:N){
cv.poly$res=sapply(cv.poly$C,function(C) cv.auc(xtrain3,ytrain3,C,type="C-svc",kernel=polydot(degree=2),scaled=c()))
Cmax=cv.poly$C[which.max(cv.poly$res)]
print(Cmax)
CC=CC+Cmax
print(i)
}
CC=CC/N

##Meilleur C moyen trouvé =  87.45796

BestPolyKernel2 = ksvm(datat,ytrain,type="C-svc",C=87.45796,kernel=polydot(degree=2),scale=c())
PolyPred2=predict(BestPolyKernel,FinalTest)







 
