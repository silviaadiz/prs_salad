library(optparse)



option_list = list(
    make_option(c("--wd"), type="character", default=NULL,help="Working directory",metavar="character"),
    make_option(c("--bd"), type="character", default=NULL,help="Binaries directory",metavar="character"),
    make_option(c("--vcf"), type="logical", default=FALSE,help="VCF file? Default BINARY",metavar="character"),
    make_option(c("--bin"), type="character", default="chr",help="Binary/vcf files pattern (in case there are several chr files in wd)",metavar="character"),
    make_option(c("--chr"), type="character", default=NULL,help="Which chr to use (separated by commas). Default 1:23",metavar="character"),
    make_option(c("--summary"), type="character", default=NULL,help="Summary statistics file",metavar="character"),
    make_option(c("--p1"), type="double", default=1e-05,help="P-val for clumping. Default 1e-05"),
    make_option(c("--p2"), type="double", default=0.01,help="Secondary threshold for clumping snps. Default 0.01"),
    make_option(c("--rsq"), type="double", default=0.5,help="LD-threshold for clumping. Default 0.5"),
    make_option(c("--kb"), type="integer", default=250,help="Kb distance for clumpin. Default 250"),
    make_option(c("--beta"), type="character", default="beta",help="Beta field on sumstats. Default beta "),
    make_option(c("--pval"), type="character", default="pval",help="Pval field on sumstats. Default pval "),
    make_option(c("--snp"), type="character", default="SNP",help="SNP field on sumstats. Default SNP "),
    make_option(c("--ea"), type="character", default="A1",help="Effect allele on sumstats. Default A1"),
    make_option(c("--pheno"), type="character", default=NULL,help="Pheno file (txt SEPARATED BY SPACES)"),
    make_option(c("--pheno_field"), type="character", default=NULL,help="Pheno field on pheno file"),
    make_option(c("--logistic"), type="logical", default=TRUE,help="Pheno is binary or continous?"),
    make_option(c("--covar"), type="character", default=NULL,help="Covariates for the model separated by comma"),
	  make_option(c("--kfold"), type="logical", default=FALSE,help="K-folds cross validation? Output is just performance measures"),
	  make_option(c("--kn"), type="integer", default=NULL,help="If K-Folds xval is true, number of K-folds"),
    make_option(c("--save"), type="logical", default=FALSE,help="If not doing cross-val, save models as RData?"),
    make_option(c("--out"), type="character", default="out",help="Output prefix for all files (analysis name)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$wd)|is.null(opt$summary)|is.null(opt$pheno)|is.null(opt$pheno_field)|is.null(opt$out)) {
  print_help(opt_parser)
  stop("Missing arguments", call.=FALSE)
}



if (opt$kfold==TRUE & is.null(opt$kn)) {
  print_help(opt_parser)
  stop("Selected K-fold cross validation but no number of K-folds", call.=FALSE)
}

pheno<-read.table(opt$pheno)
pheno["fenotipo"]<-pheno[opt$pheno_field]library(optparse)

# SILVIA ES QUE VER O DE CARGAR PLINk!

option_list = list(
    make_option(c("--wd"), type="character", default=NULL,help="Working directory",metavar="character"),
    make_option(c("--bd"), type="character", default=NULL,help="Binaries directory",metavar="character"),
    make_option(c("--vcf"), type="logical", default=FALSE,help="VCF file? Default BINARY",metavar="character"),
    make_option(c("--bin"), type="character", default="chr",help="Binary/vcf files pattern (in case there are several chr files in wd)",metavar="character"),
    make_option(c("--chr"), type="character", default=NULL,help="Which chr to use (separated by commas). Default 1:23",metavar="character"),
    make_option(c("--summary"), type="character", default=NULL,help="Summary statistics file",metavar="character"),
    make_option(c("--p1"), type="double", default=1e-05,help="P-val for clumping. Default 1e-05"),
    make_option(c("--p2"), type="double", default=0.01,help="Secondary threshold for clumping snps. Default 0.01"),
    make_option(c("--rsq"), type="double", default=0.5,help="LD-threshold for clumping. Default 0.5"),
    make_option(c("--kb"), type="integer", default=250,help="Kb distance for clumpin. Default 250"),
    make_option(c("--beta"), type="character", default="beta",help="Beta field on sumstats. Default beta "),
    make_option(c("--pval"), type="character", default="pval",help="Pval field on sumstats. Default pval "),
    make_option(c("--snp"), type="character", default="SNP",help="SNP field on sumstats. Default SNP "),
    make_option(c("--ea"), type="character", default="A1",help="Effect allele on sumstats. Default A1"),
    make_option(c("--pheno"), type="character", default=NULL,help="Pheno file (txt SEPARATED BY SPACES)"),
    make_option(c("--pheno_field"), type="character", default=NULL,help="Pheno field on pheno file"),
    make_option(c("--logistic"), type="logical", default=TRUE,help="Pheno is binary or continous?"),
    make_option(c("--covar"), type="character", default=NULL,help="Covariates for the model separated by comma"),
	  make_option(c("--kfold"), type="logical", default=FALSE,help="K-folds cross validation? Output is just performance measures"),
	  make_option(c("--kn"), type="integer", default=NULL,help="If K-Folds xval is true, number of K-folds"),
    make_option(c("--save"), type="logical", default=FALSE,help="If not doing cross-val, save models as RData?"),
    make_option(c("--out"), type="character", default="out",help="Output prefix for all files (analysis name)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$wd)|is.null(opt$summary)|is.null(opt$pheno)|is.null(opt$pheno_field)|is.null(opt$out)) {
  print_help(opt_parser)
  stop("Missing arguments", call.=FALSE)
}



if (opt$kfold==TRUE & is.null(opt$kn)) {
  print_help(opt_parser)
  stop("Selected K-fold cross validation but no number of K-folds", call.=FALSE)
}

pheno<-read.table(opt$pheno)
pheno["fenotipo"]<-pheno[opt$pheno_field]
fenotipo<-as.factor(pheno$fenotipo)

if (opt$logistic==TRUE){
		if (any(levels(fenotipo)==c("0","1"))){} else {stop("ERROR: Logistic regression is selected but column is not coded as 1/0")}}

if (opt$kfold==TRUE & opt$save==TRUE) {
  print("There is no model saving with kfold cross-val", call.=FALSE)
}

#system2("module load plink")
wd<-opt$wd

if (opt$vcf==TRUE){
list_of_files<-gsub(".vcf|.vcf.gz","",list.files(path=opt$bd,pattern=paste0("*.*",opt$bin,"*.*vcf"),full.names=T))
midf <- data.frame(archivo = list_of_files,
                   chr = NA)
                   
      for (i in 1:23){this <- which(grepl(paste0("chr",i),midf$archivo))
      midf[this,"chr"] <- i}
      
setwd(wd)
if (is.null(opt$chr)){
	for (i in seq_along(1:23)){
	bin<-midf[midf$chr==i,"archivo"]
	system2("plink",args = c(paste0("--vcf ",bin,sep="")," --allow-no-sex",
	paste0("--clump ", opt$summary,sep=""),paste0("--clump-field ", opt$pval,sep=""),
	paste0("--clump-snp-field ", opt$snp,sep=""), paste0("--clump-p1 ", opt$p1,sep=""),paste0("--clump-p2 ",opt$p2,sep=""),paste0("--clump-r2 ",opt$rsq,sep=""),paste0("--clump-kb ",opt$kb,sep=""),paste0("--out CLUMPED_",opt$out,"_",opt$p1,"_chr",i,sep="")," --silent"))
	}
 } else {
   seq_chr<-as.numeric(strsplit((opt$chr),",")[[1]])
   seq_chr<-as.vector(seq_chr)
   for (i in seq_chr){
   bin<-list_of_files[grepl(paste0("(?<![0-9])",i,"(?![0-9])"),list_of_files,perl=T)]
	system2("plink",args = c(paste0("--bfile ",bin,sep="")," --allow-no-sex",
	paste0("--clump ", opt$summary,sep=""),paste0("--clump-field ", opt$pval,sep=""),
	paste0("--clump-snp-field ", opt$snp,sep=""), paste0("--clump-p1 ", opt$p1,sep=""),paste0("--clump-p2 ",opt$p2,sep=""),paste0("--clump-r2 ",opt$rsq,sep=""),paste0("--clump-kb ",opt$kb,sep=""),paste0("--out CLUMPED_",opt$out,"_",opt$p1,"_chr",i,sep="")," --silent"))
	}
 }
} else {

list_of_files<-gsub(".bim","",list.files(path=opt$bd,pattern=paste0("*.*",opt$bin,"*.*bim"),full.names=T))
midf <- data.frame(archivo = list_of_files,
                   chr = NA)
                   
      for (i in 1:23){this <- which(grepl(paste0("chr",i),midf$archivo))
      midf[this,"chr"] <- i}

setwd(wd)
if (is.null(opt$chr)){
	for (i in 1:23){
	bin<-midf[midf$chr==i,"archivo"]
	system2("plink",args = c("--bfile ",bin," --allow-no-sex",
	paste0("--clump ", opt$summary,sep=""),paste0("--clump-field ", opt$pval,sep=""),
	paste0("--clump-snp-field ", opt$snp,sep=""), paste0("--clump-p1 ", opt$p1,sep=""),paste0("--clump-p2 ",opt$p2,sep=""),paste0("--clump-r2 ",opt$rsq,sep=""),paste0("--clump-kb ",opt$kb,sep=""),paste0("--out CLUMPED_",opt$out,"_",opt$p1,"_chr",i,sep="")," --silent"))
	}
 } else {
   seq_chr<-as.numeric(strsplit((opt$chr),",")[[1]])
   seq_chr<-as.vector(seq_chr)
   for (i in seq_chr){
   bin<-midf[midf$chr==i,"archivo"]
	system2("plink",args = c(paste0("--bfile ",bin,sep="")," --allow-no-sex",
	paste0("--clump ", opt$summary,sep=""),paste0("--clump-field ", opt$pval,sep=""),
	paste0("--clump-snp-field ", opt$snp,sep=""), paste0("--clump-p1 ", opt$p1,sep=""),paste0("--clump-p2 ",opt$p2,sep=""),paste0("--clump-r2 ",opt$rsq,sep=""),paste0("--clump-kb ",opt$kb,sep=""),paste0("--out CLUMPED_",opt$out,"_",opt$p1,"_chr",i,sep="")," --silent"))
	}
 }
}
 

 
if (is.null(list.files(pattern="CLUMPED*.clumped"))){
stop("ERROR: Check sumstats colnames or SNP names")}

list_clumped<-list.files(pattern=paste0("CLUMPED_",opt$out,"_",opt$p1,"_chr","*.*clumped"))
n=1

for (i in list_clumped){
  skip<-FALSE
  tryCatch({
  n=n+1
  clump_chr<-read.delim(i,sep="",dec=".",header=T)
  assign(paste("clumped",n,sep=""),data.frame("SNP"=clump_chr$SNP))}, error=function(e){skip<<-TRUE})
  if(skip){next}
  
}
clump_all<-do.call(rbind,mget(ls(pattern="clumped")))

write.table(clump_all,paste0("index_snps_",opt$out,"_",opt$p1,".txt",sep=""),col.names = F,row.names=F,append=F,quote=F)

rmv<-rmv<-c(paste0("CLUMPED_",opt$out,"_",opt$p1,"|(_chr)|",".clumped"))
list_extract<-gsub(rmv,"",list_clumped)

for (i in list_extract){
bin<-midf[midf$chr==i,"archivo"]
system2("plink",args = c(paste0("--bfile ",bin,sep=""),"--allow-no-sex --extract ",
paste0("index_snps_",opt$out,"_",opt$p1,".txt")," --make-bed ",
paste0("--out extracted_snps_index_",opt$out,"_",opt$p1,"_chr",i,sep="")," --silent"))}

if (is.null(list.files(pattern="extracted*.bim"))){stop("ERROR: No SNPs extracted")}

system(paste0("ls extracted_snps_index_",opt$out,"_",opt$p1,"_chr*.bim | sed \"s/.bim//g\" > merge_list",opt$out,opt$p1,".txt",sep=""))

system2("plink",args = c("--allow-no-sex ",paste0("--merge-list merge_list",opt$out,opt$p1,".txt"),
"--indiv-sort 0 ","--make-bed ",paste0("--out ALLclumped_",opt$out,"_",opt$p1,sep="")," --silent")) 



gwas_all<-read.table(opt$summary,header=T)

clumpeados_all<-read.table(paste0("index_snps_",opt$out,"_",opt$p1,".txt",sep=""))
todo_all<-merge(gwas_all,clumpeados_all,by.y="V1",by.x=opt$snp)



todo_all2<-todo_all[,c(opt$snp,opt$ea,opt$beta)] 

write.table(todo_all2,paste0("input_",opt$out,"_",opt$p1,"score.txt",sep=""),quote=F,row.names=F)

system2("plink",args = c(paste0("--bfile ALLclumped_",opt$out,"_",opt$p1,sep="")," --allow-no-sex"," --score",
paste0(" input_",opt$out,"_",opt$p1,"score.txt",sep=""), " 1 2 3",paste0(" --out PRS_",opt$out,"_",opt$p1,sep="")," --silent"))


#------ Models
prs<-read.table(paste0("PRS_",opt$out,"_",opt$p1,".profile",sep=""),header=T)

pheno2<-merge(pheno,prs,by="IID")
pheno2$st.score<-(pheno2$SCORE-mean(pheno2$SCORE))/sd(pheno2$SCORE)


library(pROC)
aucval<-as.numeric()
#sensi<-as.numeric()
#spec<-as.numeric()
aucval_prs<-as.numeric()
#sensi_prs<-as.numeric()
#spec_prs<-as.numeric()
rsq<-as.numeric()
rsq_prs<-as.numeric()
rmse<-as.numeric()
rmse_prs<-as.numeric()
k<-as.numeric()
RMSE <- function(x) { sqrt(mean(x^2)) }


set.seed(394855)

results <- list()
formula<-gsub(",","+",paste("fenotipo","~",opt$covar,sep=""))
formula2<-gsub(",","+",paste("fenotipo","~",opt$covar,"+st.score",sep=""))  
		
if (opt$kfold==TRUE){
	phenor <- pheno2[sample(nrow(pheno2)),]
	fold  <- cut(seq(1,nrow(phenor)),breaks=opt$kn,labels=FALSE)
	if (opt$logistic==TRUE){
			for(i in 1:10){
				test <- phenor[which(fold==i,arr.ind=TRUE),]
				train<-phenor[-(which(fold==i,arr.ind=TRUE)),]
				mod1<-glm(formula,train,family=binomial(link="logit"))
				preds<-predict(mod1,test,type="response")
				roc_a<-roc(test$fenotipo,predictor=preds)
				aucval[i]<-roc_a$auc
				#spec[i]<-roc_a$specificities
				#sens[i]<-roc_a$sensitivities
				
				mod2<-glm(formula2,train,family=binomial(link="logit"))
				preds2<-predict(mod2,test,type="response")
				roc_a2<-roc(test$fenotipo,predictor=preds2)
				aucval_prs[i]<-roc_a2$auc
				#spec_prs[i]<-roc_a2$specificities
				#sens_prs[i]<-roc_a2$sensitivities
				k[i]<-i

				}
				results<-cbind(k,aucval,aucval_prs)
				write.table(results,paste0(opt$out,"_LOGISTIC_",opt$p1,"_kfold",opt$kn,"_xval_performance.txt",sep=""),quote=F,row.names=F)

	}else{
			for(i in 1:10){
				test <- phenor[which(fold==i,arr.ind=TRUE),]
				train<-phenor[-(which(fold==i,arr.ind=TRUE)),]
				mod1<-lm(formula,train)
				preds<-predict(mod1,test,type="response")
				rsq[i]<-summary(mod1)$adj.r.squared
				rmse[i]<-RMSE(mod1$residuals)
				
				mod2<-lm(formula2,train)
				preds2<-predict(mod2,test,type="response")
				rsq_prs[i]<-summary(mod2)$adj.r.squared
				rmse_prs[i]<-RMSE(mod2$residuals)
				k[i]<-i
				}
				results<-cbind(k,rsq,rmse,rsq_prs,rmse_prs)
		write.table(results,paste0(opt$out,"_LINEAR_",opt$p1,"_kfold",opt$kn,"_xval_performance.txt",sep=""),quote=F,row.names=F)
		}

}else {
	n <- floor(0.70 * nrow(pheno2))
	train_ind <- sample(seq_len(nrow(pheno2)), size = n)
	train <- pheno2[train_ind, ]
	test <- pheno2[-train_ind, ]
	if (opt$logistic==TRUE){
			
			mod1<-glm(formula,train,family=binomial(link="logit"))
			preds<-predict(mod1,test,type="response")
			roc_a<-roc(test$fenotipo,predictor=preds)
			aucval<-roc_a$auc
			#spec<-roc_a$specificities
			#sens<-roc_a$sensitivities
		
			mod2<-glm(formula2,train,family=binomial(link="logit"))
			preds2<-predict(mod2,test,type="response")
			roc_a2<-roc(test$fenotipo,predictor=preds2)
			aucval_prs<-roc_a2$auc
			#spec_prs<-roc_a2$specificities
			#sens_prs<-roc_a2$sensitivities
		
		results<-cbind(aucval,aucval_prs)
		write.table(results,paste0(opt$out,"_LOGISTIC_",opt$p1,"_model_performance.txt",sep=""),quote=F,row.names=F)
		if(opt$save==TRUE){saving<-list(mod1,preds,roc_a,mod2,preds2,roc_a2)
    saveRDS(saving,file=paste0(opt$out,"_",opt$p1,"_logistic_models.rds",sep=""))}
		}else{
				mod1<-lm(formula,train)
				preds<-predict(mod1,test,type="response")
				rsq<-summary(mod1)$adj.r.squared
				rmse<-RMSE(mod1)
				
				mod2<-lm(formula2,train)
				preds2<-predict(mod2,test,type="response")
				rsq_prs<-summary(mod2)$adj.r.squared
				rmse_prs<-RMSE(mod2)
				
		results<-cbind(rsq,rmse,rsq_prs,rmse_prs)
		write.table(results,paste0(opt$out,"_LINEAR_",opt$p1,"_model_performance.txt",sep=""),quote=F,row.names=F)
		if(opt$save==TRUE){saving<-list(mod1,preds,rsq,rmse)
    saveRDS(saving,file=paste0(opt$out,"_",opt$p1,"_linear_models.rds",sep=""))}
		}
}
fenotipo<-as.factor(pheno$fenotipo)

if (opt$logistic==TRUE){
		if (any(levels(fenotipo)==c("0","1"))){} else {stop("ERROR: Logistic regression is selected but column is not coded as 1/0")}}

if (opt$kfold==TRUE & opt$save==TRUE) {
  print("There is no model saving with kfold cross-val", call.=FALSE)
}

#system2("module load plink")
wd<-opt$wd

if (opt$vcf==TRUE){
list_of_files<-gsub(".vcf|.vcf.gz","",list.files(path=opt$bd,pattern=paste0("*.*",opt$bin,"*.*vcf"),full.names=T))
midf <- data.frame(archivo = list_of_files,
                   chr = NA)
                   
      for (i in 1:23){this <- which(grepl(paste0("chr",i),midf$archivo))
      midf[this,"chr"] <- i}
      
setwd(wd)
if (is.null(opt$chr)){
	for (i in seq_along(1:23)){
	bin<-midf[midf$chr==i,"archivo"]
	system2("plink",args = c(paste0("--vcf ",bin,sep="")," --allow-no-sex",
	paste0("--clump ", opt$summary,sep=""),paste0("--clump-field ", opt$pval,sep=""),
	paste0("--clump-snp-field ", opt$snp,sep=""), paste0("--clump-p1 ", opt$p1,sep=""),paste0("--clump-p2 ",opt$p2,sep=""),paste0("--clump-r2 ",opt$rsq,sep=""),paste0("--clump-kb ",opt$kb,sep=""),paste0("--out CLUMPED_",opt$out,"_",opt$p1,"_chr",i,sep="")," --silent"))
	}
 } else {
   seq_chr<-as.numeric(strsplit((opt$chr),",")[[1]])
   seq_chr<-as.vector(seq_chr)
   for (i in seq_chr){
   bin<-list_of_files[grepl(paste0("(?<![0-9])",i,"(?![0-9])"),list_of_files,perl=T)]
	system2("plink",args = c(paste0("--bfile ",bin,sep="")," --allow-no-sex",
	paste0("--clump ", opt$summary,sep=""),paste0("--clump-field ", opt$pval,sep=""),
	paste0("--clump-snp-field ", opt$snp,sep=""), paste0("--clump-p1 ", opt$p1,sep=""),paste0("--clump-p2 ",opt$p2,sep=""),paste0("--clump-r2 ",opt$rsq,sep=""),paste0("--clump-kb ",opt$kb,sep=""),paste0("--out CLUMPED_",opt$out,"_",opt$p1,"_chr",i,sep="")," --silent"))
	}
 }
} else {

list_of_files<-gsub(".bim","",list.files(path=opt$bd,pattern=paste0("*.*",opt$bin,"*.*bim"),full.names=T))
midf <- data.frame(archivo = list_of_files,
                   chr = NA)
                   
      for (i in 1:23){this <- which(grepl(paste0("chr",i),midf$archivo))
      midf[this,"chr"] <- i}

setwd(wd)
if (is.null(opt$chr)){
	for (i in 1:23){
	bin<-midf[midf$chr==i,"archivo"]
	system2("plink",args = c("--bfile ",bin," --allow-no-sex",
	paste0("--clump ", opt$summary,sep=""),paste0("--clump-field ", opt$pval,sep=""),
	paste0("--clump-snp-field ", opt$snp,sep=""), paste0("--clump-p1 ", opt$p1,sep=""),paste0("--clump-p2 ",opt$p2,sep=""),paste0("--clump-r2 ",opt$rsq,sep=""),paste0("--clump-kb ",opt$kb,sep=""),paste0("--out CLUMPED_",opt$out,"_",opt$p1,"_chr",i,sep="")," --silent"))
	}
 } else {
   seq_chr<-as.numeric(strsplit((opt$chr),",")[[1]])
   seq_chr<-as.vector(seq_chr)
   for (i in seq_chr){
   bin<-midf[midf$chr==i,"archivo"]
	system2("plink",args = c(paste0("--bfile ",bin,sep="")," --allow-no-sex",
	paste0("--clump ", opt$summary,sep=""),paste0("--clump-field ", opt$pval,sep=""),
	paste0("--clump-snp-field ", opt$snp,sep=""), paste0("--clump-p1 ", opt$p1,sep=""),paste0("--clump-p2 ",opt$p2,sep=""),paste0("--clump-r2 ",opt$rsq,sep=""),paste0("--clump-kb ",opt$kb,sep=""),paste0("--out CLUMPED_",opt$out,"_",opt$p1,"_chr",i,sep="")," --silent"))
	}
 }
}
 

 
if (is.null(list.files(pattern="CLUMPED*.clumped"))){
stop("ERROR: Check sumstats colnames or SNP names")}

list_clumped<-list.files(pattern=paste0("CLUMPED_",opt$out,"_",opt$p1,"_chr","*.*clumped"))
n=1

for (i in list_clumped){
  skip<-FALSE
  tryCatch({
  n=n+1
  clump_chr<-read.delim(i,sep="",dec=".",header=T)
  assign(paste("clumped",n,sep=""),data.frame("SNP"=clump_chr$SNP))}, error=function(e){skip<<-TRUE})
  if(skip){next}
  
}
clump_all<-do.call(rbind,mget(ls(pattern="clumped")))

write.table(clump_all,paste0("index_snps_",opt$out,"_",opt$p1,".txt",sep=""),col.names = F,row.names=F,append=F,quote=F)

rmv<-rmv<-c(paste0("CLUMPED_",opt$out,"_",opt$p1,"|(_chr)|",".clumped"))
list_extract<-gsub(rmv,"",list_clumped)

for (i in list_extract){
bin<-midf[midf$chr==i,"archivo"]
system2("plink",args = c(paste0("--bfile ",bin,sep=""),"--allow-no-sex --extract ",
paste0("index_snps_",opt$out,"_",opt$p1,".txt")," --make-bed ",
paste0("--out extracted_snps_index_",opt$out,"_",opt$p1,"_chr",i,sep="")," --silent"))}

if (is.null(list.files(pattern="extracted*.bim"))){stop("ERROR: No SNPs extracted")}

system(paste0("ls extracted_snps_index_",opt$out,"_",opt$p1,"_chr*.bim | sed \"s/.bim//g\" > merge_list",opt$out,opt$p1,".txt",sep=""))

system2("plink",args = c("--allow-no-sex ",paste0("--merge-list merge_list",opt$out,opt$p1,".txt"),
"--indiv-sort 0 ","--make-bed ",paste0("--out ALLclumped_",opt$out,"_",opt$p1,sep="")," --silent")) 



gwas_all<-read.table(opt$summary,header=T)

clumpeados_all<-read.table(paste0("index_snps_",opt$out,"_",opt$p1,".txt",sep=""))
todo_all<-merge(gwas_all,clumpeados_all,by.y="V1",by.x=opt$snp)



todo_all2<-todo_all[,c(opt$snp,opt$ea,opt$beta)] 

write.table(todo_all2,paste0("input_",opt$out,"_",opt$p1,"score.txt",sep=""),quote=F,row.names=F)

system2("plink",args = c(paste0("--bfile ALLclumped_",opt$out,"_",opt$p1,sep="")," --allow-no-sex"," --score",
paste0(" input_",opt$out,"_",opt$p1,"score.txt",sep=""), " 1 2 3",paste0(" --out PRS_",opt$out,"_",opt$p1,sep="")," --silent"))


#------ Models
prs<-read.table(paste0("PRS_",opt$out,"_",opt$p1,".profile",sep=""),header=T)

pheno2<-merge(pheno,prs,by="IID")
pheno2$st.score<-(pheno2$SCORE-mean(pheno2$SCORE))/sd(pheno2$SCORE)


library(pROC)
aucval<-as.numeric()
#sensi<-as.numeric()
#spec<-as.numeric()
aucval_prs<-as.numeric()
#sensi_prs<-as.numeric()
#spec_prs<-as.numeric()
rsq<-as.numeric()
rsq_prs<-as.numeric()
rmse<-as.numeric()
rmse_prs<-as.numeric()
k<-as.numeric()
RMSE <- function(x) { sqrt(mean(x^2)) }


set.seed(394855)

results <- list()
formula<-gsub(",","+",paste("fenotipo","~",opt$covar,sep=""))
formula2<-gsub(",","+",paste("fenotipo","~",opt$covar,"+st.score",sep=""))  
		
if (opt$kfold==TRUE){
	phenor <- pheno2[sample(nrow(pheno2)),]
	fold  <- cut(seq(1,nrow(phenor)),breaks=opt$kn,labels=FALSE)
	if (opt$logistic==TRUE){
			for(i in 1:10){
				test <- phenor[which(fold==i,arr.ind=TRUE),]
				train<-phenor[-(which(fold==i,arr.ind=TRUE)),]
				mod1<-glm(formula,train,family=binomial(link="logit"))
				preds<-predict(mod1,test,type="response")
				roc_a<-roc(test$fenotipo,predictor=preds)
				aucval[i]<-roc_a$auc
				#spec[i]<-roc_a$specificities
				#sens[i]<-roc_a$sensitivities
				
				mod2<-glm(formula2,train,family=binomial(link="logit"))
				preds2<-predict(mod2,test,type="response")
				roc_a2<-roc(test$fenotipo,predictor=preds2)
				aucval_prs[i]<-roc_a2$auc
				#spec_prs[i]<-roc_a2$specificities
				#sens_prs[i]<-roc_a2$sensitivities
				k[i]<-i

				}
				results<-cbind(k,aucval,aucval_prs)
				write.table(results,paste0(opt$out,"_LOGISTIC_",opt$p1,"_kfold",opt$kn,"_xval_performance.txt",sep=""),quote=F,row.names=F)
	################### PROBLEMAS CO SAVING RDS ######################
	}else{
			for(i in 1:10){
				test <- phenor[which(fold==i,arr.ind=TRUE),]
				train<-phenor[-(which(fold==i,arr.ind=TRUE)),]
				mod1<-lm(formula,train)
				preds<-predict(mod1,test,type="response")
				rsq[i]<-summary(mod1)$adj.r.squared
				rmse[i]<-RMSE(mod1$residuals)
				
				mod2<-lm(formula2,train)
				preds2<-predict(mod2,test,type="response")
				rsq_prs[i]<-summary(mod2)$adj.r.squared
				rmse_prs[i]<-RMSE(mod2$residuals)
				k[i]<-i
				}
				results<-cbind(k,rsq,rmse,rsq_prs,rmse_prs)
		write.table(results,paste0(opt$out,"_LINEAR_",opt$p1,"_kfold",opt$kn,"_xval_performance.txt",sep=""),quote=F,row.names=F)
		}

}else {
	n <- floor(0.70 * nrow(pheno2))
	train_ind <- sample(seq_len(nrow(pheno2)), size = n)
	train <- pheno2[train_ind, ]
	test <- pheno2[-train_ind, ]
	if (opt$logistic==TRUE){
			
			mod1<-glm(formula,train,family=binomial(link="logit"))
			preds<-predict(mod1,test,type="response")
			roc_a<-roc(test$fenotipo,predictor=preds)
			aucval<-roc_a$auc
			#spec<-roc_a$specificities
			#sens<-roc_a$sensitivities
		
			mod2<-glm(formula2,train,family=binomial(link="logit"))
			preds2<-predict(mod2,test,type="response")
			roc_a2<-roc(test$fenotipo,predictor=preds2)
			aucval_prs<-roc_a2$auc
			#spec_prs<-roc_a2$specificities
			#sens_prs<-roc_a2$sensitivities
		
		results<-cbind(aucval,aucval_prs)
		write.table(results,paste0(opt$out,"_LOGISTIC_",opt$p1,"_model_performance.txt",sep=""),quote=F,row.names=F)
		if(opt$save==TRUE){saving<-list(mod1,preds,roc_a,mod2,preds2,roc_a2)
    saveRDS(saving,file=paste0(opt$out,"_",opt$p1,"_logistic_models.rds",sep=""))}
		}else{
				mod1<-lm(formula,train)
				preds<-predict(mod1,test,type="response")
				rsq<-summary(mod1)$adj.r.squared
				rmse<-RMSE(mod1)
				
				mod2<-lm(formula2,train)
				preds2<-predict(mod2,test,type="response")
				rsq_prs<-summary(mod2)$adj.r.squared
				rmse_prs<-RMSE(mod2)
				
		results<-cbind(rsq,rmse,rsq_prs,rmse_prs)
		write.table(results,paste0(opt$out,"_LINEAR_",opt$p1,"_model_performance.txt",sep=""),quote=F,row.names=F)
		if(opt$save==TRUE){saving<-list(mod1,preds,rsq,rmse)
    saveRDS(saving,file=paste0(opt$out,"_",opt$p1,"_linear_models.rds",sep=""))}
		}
}
