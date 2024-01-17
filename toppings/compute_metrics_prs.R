 library(optparse)


option_list = list(
    make_option(c("--wd"), type="character", default=NULL,help="Working directory",metavar="character"),
  make_option(c("--score"),type="character",default=NULL),
  make_option(c("--p1"), type="double", default=1e-05,help="P-val for clumping. Default 1e-05"),
  make_option(c("--rsq"), type="double", default=0.5,help="LD-threshold for clumping. Default 0.5"),
  make_option(c("--kb"), type="integer", default=250,help="Kb distance for clumpin. Default 250"),
  make_option(c("--pheno"), type="character", default=NULL,help="Pheno file (txt SEPARATED BY SPACES)"),
  make_option(c("--pheno_field"), type="character", default=NULL,help="Pheno field on pheno file"),
  make_option(c("--covar"), type="character", default=NULL,help="Covariates for the model separated by comma"),
  make_option(c("--nrep"), type="integer", default=NULL,help="Bootstrap reps"),
  make_option(c("--out"), type="character", default="out",help="Output prefix for all files (analysis name)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
wd<-opt$wd
setwd(wd)
pheno<-read.table(opt$pheno,header=T)
  pheno["fenotipo"]<-pheno[opt$pheno_field]
  pheno$fenotipo<-as.factor(pheno$fenotipo)
  #------ Models
  prs<-read.table(paste0(opt$score,"_",opt$p1,"_",opt$rsq,"_",opt$kb,".profile",sep=""),header=T)
  
  pheno2<-merge(pheno,prs,by="IID")
  pheno2$st.score<-(pheno2$SCORE-mean(pheno2$SCORE))/sd(pheno2$SCORE)
  

  library(rms)
  library(boot)
  set.seed(394855)
  reps<-opt$nrep

  f1<-as.formula(gsub(",","+",paste("fenotipo","~",opt$covar,sep="")))
  f2<-as.formula(gsub(",","+",paste("fenotipo","~",opt$covar,"+st.score",sep="")))
  
  dif_rsq <- function(formula1,formula2, data, index) {
   bt<- data[index,] 
  mod1 <-lrm(formula1,bt)
  mod2 <-lrm(formula2,bt)
  dif<-mod2$stats["R2"]-mod1$stats["R2"]

  return(dif)}
  
  results <- boot(data=pheno2, statistic=dif_rsq,
   R=reps, formula1=f1,formula2=f2)
outp<-mean(results$t)
cis<-boot.ci(results,type="basic")

r2_iclow<-cis$basic[4]
r2_ichigh<-cis$basic[5]
outp2<-cbind(outp,r2_iclow,r2_ichigh)
names(outp2)<-c("r2_inc","ic_low","ic_high")

outp2$PRS<-paste0(opt$p1,"_",opt$rsq,"_",opt$kb,sep="")

write.table(outp2,paste0(opt$out,"_",opt$p1,"_",opt$rsq,"_",opt$kb,"_","validation_bootstrap_",opt$nrep,".txt",sep=""),quote=F)

warnings()
