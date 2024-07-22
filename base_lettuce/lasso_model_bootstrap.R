library(glmnet)
library(pROC)

boot_lasso<-function(data,pheno,covar,interactions=F,interactions_covar=NULL,B){
set.seed(1711)
dat_na<-data[complete.cases(data[,c(as.character(pheno),as.character(covar))]),,drop=F]
  
if(interactions==T&!is.null(interactions_covar)){
f1<-as.formula(sprintf('%s ~ %s+%s',as.character(pheno),paste(covar,collapse="+"),paste(interactions_covar,collapse="+")))
}else if(interactions==F){
f1<-as.formula(sprintf('%s ~ %s',as.character(pheno),paste(covar,collapse="+")))}
  
x<-model.matrix(f1,data=dat_na)[,-1]
n<-nrow(dat_na)
B<-B
coef_bootstrap_1se=matrix(NA,nrow=B,ncol=ncol(x)+1)
colnames(coef_bootstrap_1se)=c("intercept",colnames(x))
boot.lambda=matrix(NA,nrow=B,ncol=1)
colnames(boot.lambda)=c("1se")
res<-list()
preds<-matrix(nrow=nrow(dat_na),ncol=B)
colnames(preds)<-c(seq(1,B,1))
aucs<-matrix(NA,nrow=B,ncol=1)
  
for(i in 1:B){
ind=sample(1:n,size=n,replace=T)#bootstrap index
x_boot<-x[ind,]
y_boot<-dat_na[ind,pheno]
boot_cv_lasso=cv.glmnet(x=x_boot,y=y_boot,alpha=1)#calculamos lambda 
coef_bootstrap_1se[i,]=as.numeric(coef(boot_cv_lasso$glmnet.fit,s=boot_cv_lasso$lambda.1se))
boot.lambda[i,]=c(boot_cv_lasso$lambda.1se)
preds_boot<-predict(boot_cv_lasso,newx=x_boot,type="response",s="lambda.1se")
#Predictions in the test sample
preds[,i]<-predict(boot_cv_lasso,newx=x,type="response",s="lambda.1se")
#Get optimism
roc_curve_boot<-auc(roc(y_boot,as.vector(preds_boot)))
roc_curve_test<-auc(roc(dat_na[,pheno],as.vector(preds[,i])))
aucs[i,]<-roc_curve_boot-roc_curve_test
}
  
prob.nonzero.1se=apply(abs(coef_bootstrap_1se)>1e-5,2,mean)
res$coef<-coef_bootstrap_1se;res$lambda<-boot.lambda;res$prob.nonzero<-prob.nonzero.1se;res$pred_cal<-preds;res$opt_cor_auc<-aucs
return(res)
}

