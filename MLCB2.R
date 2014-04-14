#téléchargement des données

xtrain = read.table("//paradis/eleves/CDell.Aiera/Bureau/MLCB/xtrain.txt")
ytrain = read.table("//paradis/eleves/CDell.Aiera/Bureau/MLCB/ytrain.txt")

matrix(xtrain)
data=as.matrix(xtrain[,2:185])
datat=t(data)
xtrain2=as.data.frame(datat)

#méthodes de classification par différents noyaux
install.packages("kernlab")
library(kernlab)

svp <- ksvm(datat,ytrain,type="C-svc",kernel="rbfdot",kpar=list(sigma=1),scale=FALSE)
plot(svp,data=as.data.frame(datat))
title("SVM classification plot")
#predict(object, newdata, type = "response", coupler = "minpair")


svp <- ksvm(xtrain2,ytrain,type="C-svc",kernel="vanilladot")
plot(svp,data=datat)

svp <- ksvm(T,ytrain,type="C-svc",kernel=polydot(degree=2))
plot(svp,data=x)

#summary(xtrain)