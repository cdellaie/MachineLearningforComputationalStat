library(kernlab)
library(ROCR)
library(igraph)

#téléchargement des données

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
#essai sur les set d'entraînements avec plusieurs noyaux
#le meilleur sera choisi pour être entrainé sur la base complète et prévoir sur les données à traiter
#tanhdot,laplacedot,besseldot,anovadot
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

#############################
## Fonctions préliminaires ##
#############################

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

# We compute the prediction vector by cross-validation

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



 
