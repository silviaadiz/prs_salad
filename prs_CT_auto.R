library(optparse)

option_list = list(
    make_option(c("-d","--dir"), type="character", default=NULL,help="Working directory",metavar="character"),
    make_option(c("-bin","--binary"), type="character", default="chr",help="Binary files pattern (in case there are several chr files in wd)",metavar="character"),
    make_option(c("-s","--summary"), type="character", default=NULL,help="Summary statistics file",metavar="character"),
    make_option(c("-p1","--p1"), type="double", default=1e-05,help="P-val for clumping"),
    make_option(c("-p2","--p2"), type="double", default=0.01,help="Secondary threshold for clumping snps"),
    make_option(c("-r2","--rsq"), type="double", default=0.5,help="LD-threshold for clumping"),
    make_option(c("-kb","--clump_kb"), type="integer", default=250,help="Kb distance for clumping"),
    make_option(c("-beta","--beta"), type="character", default="beta",help="Beta field on sumstats"),
    make_option(c("-pval","--pval"), type="character", default="pval",help="Pval field on sumstats"),
    make_option(c("-snp","--snp"), type="character", default="SNP",help="SNP field on sumstats"),
    make_option(c("-ea","--effect_allele"), type="character", default="A1",help="Effect allele on sumstats"),
    make_option(c("-pheno","--pheno"), type="character", default=NULL,help="Pheno file"),
    make_option(c("-pf","--pheno_field"), type="character", default=NULL,help="Pheno field on pheno file"),
    make_option(c("-logistic","--logistic"), type="logical", default=TRUE,help="Pheno is binary or continous?"),
    make_option(c("-covar","--covar"), type="character", default=NULL,help="Covariates for the model"),
    make_option(c("-save","--save_model"), type="logical", default=FALSE,help="Save models as RData?"),
    make_option(c("-out","--out"), type="character", default="out",help="Output prefix for all files (analysis name)"),
    make_option(c("-outd","--outd"), type="character", default=NULL,help="Output directory for jpegs and models, optional"),
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$dir)|is.null(opt$summary)|is.null(opt$pheno)|is.null(opt$pheno_field)|is.null(opt$out)) {
  print_help(opt_parser)
  stop("Missing arguments", call.=FALSE)
}



system2("module load plink")
wd<-opt$dir
setwd(wd)

list_of_files<-gsub(".bim","",list.files(pattern=paste0(pat,"*.bim")))


	for (i in 1:23){
 bin<-list_of_files[i]
	system2("plink",args = c(paste0("--bfile ",bin,sep=""),"--allow-no-sex",
	paste0("--clump ", opt$summary,sep=""),paste0("--clump-field ", opt$pval,sep=""),
	paste0("--clump-snp-field ", opt$snp,sep=""), paste0("--clump-p1 ", opt$p1,sep=""),paste0("--clump-p2",opt$p2,sep=""),paste0("--clump-r2",opt$rsq,sep=""),paste0("--clump-kb",opt$clump_kb,sep=""),paste0("--out CLUMPED_",opt$out,"_",opt$p1,"_chr",i,sep=""),"--silent"))
	}
 
if (is.null(list.files(pattern="CLUMPED*.clumped"))){
stop("ERROR: Check sumstats colnames or SNP names",.call=F)}

for (i in 1:23){
  skip<-FALSE
  tryCatch({
  clump_chr<-read.delim(paste("CLUMPED_",opt$out,"_",opt$p1,"_chr",i,".clumped",sep=""),sep="",dec=".",header=T)
  assign(paste("clumped",i,sep=""),data.frame("SNP"=clump_chr$SNP))}, error=function(e){skip<<-TRUE})
  if(skip){next}
  
}
clump_all<-do.call(rbind,mget(ls(pattern="clumped")))

write.table(clump_all,paste0("index_snps_",opt$out,"_",opt$p1,".txt",sep=""),col.names = F,row.names=F,append=F,quote=F)

for (i in 1:23){
bin<-list_of_files[i]
system2("plink",args = c(paste0("--bfile ",bin,sep=""),"--allow-no-sex","--extract ",
paste0("index_snps_",opt$out,"_",opt$p1,".txt"),"--make-bed",
paste0("--out extracted_snps_index_",opt$out,"_",opt$p1,"_chr",i,sep=""),"--silent"))}

if (is.null(list.files(pattern="extracted*.bim"))){stop("ERROR: No SNPs extracted",.call=F)}

system(paste0("ls extracted_snps_index_",opt$out,"_",opt$p1,"_chr*.bim | sed \"s/.bim//g\" > merge_list",opt$out,opt$p1,".txt",sep=""))

system2("plink",args = c("--allow-no-sex",paste0("--merge-list merge_list",opt$out,opt$p1,".txt"),
"--indiv-sort 0","--make-bed",paste0("--out ALLclumped_",opt$out,"_",opt$p1,sep=""),"--silent")) 



gwas_all<-read.delim(opt$summary,header=T)

clumpeados_all<-read.table(paste0("index_snps_",opt$out,"_",opt$p1,".txt",sep=""))
todo_all<-merge(gwas_all,clumpeados_all,by.y="V1",by.x=opt$snp)



todo_all2<-todo_all[,c(opt$snp,opt$effect_allele,opt$beta)] 

write.table(todo_all2,paste0("input_",opt$out,"_",opt$p1,"score.txt",sep=""),quote=F,row.names=F)

system2("plink",args = c(paste0("--bfile ALLclumped_",opt$out,"_",opt$p1,sep=""),"--allow-no-sex","--score",
paste0("input_",opt$out,"_",opt$p1,"score.txt",sep=""), "1 2 3",paste0("--out PRS_",opt$out,"_",opt$p1,sep="")))


#------ Models

if (opt$