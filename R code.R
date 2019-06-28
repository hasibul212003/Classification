data<- read.csv("data.csv")
labels <- read.csv("labels.csv")
library(plyr)
fulldata<- join(data,labels)
fulldata<- fulldata[,-1]

#Creating labels and adding it to the dataframe
names(sort(table(fulldata[,20532])))
class1indices<- grep("KIRC|BRCA",fulldata[,20532])
class<- rep(0,nrow(fulldata))
class[class1indices]<- 1
class<- as.factor(class)
fulldata<- data.frame(fulldata,class)
fulldata<- fulldata[,-20532]
#Split the data into training and test prior to cross validation/feature selection
set.seed(32603)
train<- sample(1:nrow(fulldata),0.8*nrow(fulldata),FALSE)
trainfulldata<- fulldata[train,]
testfulldata<- fulldata[-train,]




#Cross Validation for Lasso
set.seed(32603)
library(glmnet)
crossval_lasso<- cv.glmnet(x=as.matrix(trainfulldata[,1:20531]),y=trainfulldata[,20532],alpha=1,type.measure = "class",family="binomial")
optimumlambda<- crossval_lasso$lambda.min
optlambdaindex<- which(crossval_lasso$lambda==crossval_lasso$lambda.min)
plot(crossval_lasso$lambda,crossval_lasso$cvm,type="l",main="crossvalidation error as a function of lambda",xlab="lambda",ylab="cross validation error")
points(crossval_lasso$lambda[optlambdaindex],crossval_lasso$cvm[optlambdaindex],cex=3.5,col="red")


#Misclassification error for optimum lambda on the test dataset
set.seed(32603)
train_lasso<- glmnet(x=as.matrix(trainfulldata[,1:20531]),y=trainfulldata[,20532],alpha=1,lambda=optimumlambda,family="binomial")
pred_test<- predict(train_lasso,newx = as.matrix(testfulldata[,1:20531]),type="response")
pred_class<-rep(0,nrow(testfulldata))
pred_class[pred_test>0.5]<-1  
test_error<- 1-mean(pred_class==testfulldata[,20532])

#Fitting the Lasso model on the entire dataset
set.seed(32603)
y_factor<- fulldata[,20532]
lasso_mod<- glmnet(as.matrix(fulldata[,1:20531]),y=y_factor,alpha=1,family="binomial",lambda = optimumlambda)

# The non zero coefficients
coefvector<-coef(lasso_mod)[which(coef(lasso_mod) != 0)]
nonzerocoefindex<- which(coef(lasso_mod) != 0)
names(coefvector)<-c("intercept",colnames(fulldata)[nonzerocoefindex[-1]])



#Feature selection using Geneselector and Limma
library(GeneSelector)  #ran on Hiper Gator
library(limma)
x<-t(trainfulldata[,-20532])
y<-trainfulldata[,20532]
limma <- RankingLimma(x, y, type="unpaired")
top500genes<-get(load("top500genes.Rdata"))

#Split the dataset into testing and training on the selected 500 genes
selecttraindata<- data.frame(trainfulldata[,top500genes$index],trainfulldata[,20532])
selecttestdata<- data.frame(testfulldata[,top500genes$index],testfulldata[,20532])
names(selecttraindata)[501]<-"class"
names(selecttestdata)[501]<-"class"
testclass<- selecttestdata[,501]




#Testing whether variance is different in different in the two classes
variancevec<-rep(0,500)
for(i in 1:500){
  variancevec[i]<-bartlett.test(selecttraindata[,i],selecttraindata[,501])$p.value
}
#Performing qda since the pvalue for the test of equality of variance between the two groups was rejected for many genes 
library(MASS)
mylda<- qda(class~.,data=selecttraindata)
ldaclass<- predict(mylda,newdata = selecttestdata)
ldatesterror<- 1-mean(ldaclass$class==testclass)






#Plot of top two genes with response in colors
index<- top500genes$index
library(ggplot2)
mydata<-data.frame(fulldata[,index[1]],fulldata[,index[2]],fulldata[,20532])
names(mydata)<- c("gene17317","gene8032","class")
ggplot(mydata)+aes(x=mydata$gene17317,y=mydata$gene8032,color=mydata$class)+geom_point()+labs(title = "Top two genes plotted with response shown in color",x="gene17317",y="gene8032")


#svm with linear kernel
#cross validation on the training data of selected 500 genes
svm_cvlinear<- tune(svm,class~.,data=selecttraindata,kernel="linear",ranges = list(cost=c(0.001,0.01,0.1,1,5,10,100)))
bestsvmlinear<- svm_cvlinear$best.model
#prediction on test data
svmlinearpred<- predict(bestsvmlinear,newdata=selecttestdata)
svmlinearerror<- 1-mean(svmlinearpred==testclass)
