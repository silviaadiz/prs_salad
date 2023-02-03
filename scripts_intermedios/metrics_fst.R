library(optparse)


option_list = list(
  make_option(c("--binfile"),type="character",default=NULL,help="Binary LD-pruned file"),
   make_option(c("--prs"),type="character",default=NULL,help="Prs input file"),
  make_option(c("--p1"), type="double", default=1e-05,help="P-val for clumping. Default 1e-05"),
  make_option(c("--rsq"), type="double", default=0.5,help="LD-threshold for clumping. Default 0.5"),
  make_option(c("--kb"), type="integer", default=250,help="Kb distance for clumpin. Default 250"),
  make_option(c("--pheno"), type="character", default=NULL,help="Pheno file (txt SEPARATED BY SPACES)"),
  make_option(c("--pheno_field"), type="character", default=NULL,help="Pheno field on pheno file"),
  make_option(c("--logistic"), type="logical", default=TRUE,help="Pheno is binary or continous?"),
  make_option(c("--covar"), type="character", default=NULL,help="Covariates for the model separated by comma"),
  make_option(c("--kn"), type="integer", default=NULL,help="If K-Folds xval is true, number of K-folds"),
  make_option(c("--rept"), type="integer", default=NULL,help="If K-Folds xval is true, number of K-folds repetition"),
  make_option(c("--out"), type="character", default="out",help="Output prefix for all files (analysis name)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

pheno<-read.table(opt$pheno,header=T)
  pheno["fenotipo"]<-pheno[opt$pheno_field]
  pheno$fenotipo<-as.factor(pheno$fenotipo)
  #------ Models
  prs<-read.table(paste0(opt$prs,"_",opt$p1,"_",opt$rsq,"_",opt$kb,".profile",sep=""),header=T)
  
  pheno2<-merge(pheno,prs,by="IID")
  pheno2$st.score<-(pheno2$SCORE-mean(pheno2$SCORE))/sd(pheno2$SCORE)
  
  cat("ben1")
  library(rms)
  library(pROC)
  
  
  cv<-function(data,kn,r,logistic){
  r2<-as.numeric()
  r2_prs<-as.numeric()
  rsq<-as.numeric()
  rsq_prs<-as.numeric()
  k<-as.numeric()
  fst_trt_k<-as.numeric()
  fst_caco_k<-as.numeric()
 
  res_rep<-NULL
  
  rept=1
    while(rept<=r){
  datar <- data[sample(nrow(data)),] # shuffling
    fold  <- cut(seq(1,nrow(datar)),breaks=kn,labels=FALSE)
      for(i in 1:kn){
      
      # TRAIN-TEST division
        test <- datar[which(fold==i,arr.ind=TRUE),]
        train<-datar[-(which(fold==i,arr.ind=TRUE)),]
      # FST CALCULATIONS: we need to create a cluster file for plink.
      test$fst_cluster<-paste0("1")
      train$fst_cluster<-paste0("2")
      fst_export_a<-rbind(test[c("IID","IID","fst_cluster")],train[c("IID","IID","fst_cluster")])
      fst_export_b<-train[c("IID","IID","fenotipo")]
      write.table(fst_export_a,"clusters_fst_trt.txt",quote=F,row.names=F,col.names=F)
       write.table(fst_export_b,"clusters_fst_caco.txt",quote=F,row.names=F,col.names=F)
         system2("/mnt/lustre/scratch/nlsas/home/usc/gb/sdd/plink2",args=c(paste0("--bfile ", opt$binfile, " --allow-no-sex --fst CATPHENO --within clusters_fst_trt.txt --out fst_trt --silent")))
       system2("/mnt/lustre/scratch/nlsas/home/usc/gb/sdd/plink2",args=c(paste0("--bfile ", opt$binfile, " --allow-no-sex --fst CATPHENO --within clusters_fst_caco.txt --out fst_caco --silent")))
       # We use the Hudson method
       fst_trt<-read.table("fst_trt.fst.summary",header=F)
       fst_trt_k[i]<-fst_trt["V3"]
       fst_caco<-read.table("fst_caco.fst.summary",header=F)
       fst_caco_k[i]<-fst_caco["V3"]
       
       
        if(logistic==TRUE){
        mod1<-lrm(formula1,train)
        r2[i]<-mod1$stats["R2"]


               
        mod2<-lrm(formula2,train)
        r2_prs[i]<-mod2$stats["R2"]
        k[i]<-i
        
      } else {mod1<-lm(formula1,train)
        preds<-predict(mod1,test,type="response")
        rsq[i]<-summary(mod1)$adj.r.squared
        
        mod2<-lm(formula2,train)
        preds2<-predict(mod2,test,type="response")
        rsq_prs[i]<-summary(mod2)$adj.r.squared
        k[i]<-i
        
      }
      }
      if(logistic==TRUE) {results_kfold<-cbind(k,r2,r2_prs,fst_caco_k,fst_trt_k)} else {results_kfold<-cbind(k,rsq,rsq_prs,fst_caco_k,fst_trt_k)}
      res_rep<-rbind(res_rep,results_kfold)
      rept=rept+1}
      return(res_rep)
      print(res_rep)}
      
      
  set.seed(394855)
  

  formula1<-as.formula(gsub(",","+",paste("fenotipo","~",opt$covar,sep="")))
  formula2<-as.formula(gsub(",","+",paste("fenotipo","~",opt$covar,"+st.score",sep="")))

  
  k=opt$kn

res_kf<-cv(pheno2,kn=k,r=opt$rept,logistic=opt$logistic)


write.table(res_kf,paste0(opt$out,"_",opt$p1,"_",opt$rsq,"_",opt$kb,"_kfold",opt$kn,"fst.txt",sep=""),quote=F,row.names=F)

  
#}
warnings()
