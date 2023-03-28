 library(optparse)



option_list = list(
  make_option(c("--wd"), type="character", default=NULL,help="Working directory",metavar="character"),
  make_option(c("--bd"), type="character", default=NULL,help="Binaries directory",metavar="character"),
  make_option(c("--bin"), type="character", default="chr",help="Binary/vcf files pattern (in case there are several chr files in wd)",metavar="character"),
  make_option(c("--summary"), type="character", default=NULL,help="Summary statistics file",metavar="character"),
  make_option(c("--p1"), type="double", default=1e-05,help="P-val for clumping. Default 1e-05"),
  make_option(c("--p2"), type="double", default=0.01,help="P-val 2 for clumping. Default 0.01"),
  make_option(c("--rsq"), type="double", default=0.1,help="LD-threshold for clumping. Default 0.5"),
  make_option(c("--kb"), type="integer", default=250,help="Kb distance for clumpin. Default 250"),
  make_option(c("--pval"), type="character", default="pval",help="Pval field on sumstats. Default pval "),
  make_option(c("--snp"), type="character", default="SNP",help="SNP field on sumstats. Default SNP "),
  make_option(c("--out"), type="character", default="out",help="Output prefix for all files (analysis name)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

 pk<-"module load plink"
system(pk)
wd<-opt$wd
 
 setwd(wd)

  list_of_files<-gsub(".bim","",list.files(path=opt$bd,pattern=paste0(".*",opt$bin,"*.*bim"),full.names=T))
  midf <- data.frame(archivo = list_of_files,
                     chr = NA)
  
  for (i in 1:23){this <- which(grepl(paste0("chr",i),midf$archivo))
  midf[this,"chr"] <- i}

    for (i in 1:23){
      bin<-midf[midf$chr==i,"archivo"]
      system2("plink",args = c(paste0("--bfile ",bin," --allow-no-sex"),
                               paste0("--clump ", opt$summary,sep=""),paste0("--clump-field ", opt$pval,sep=""),
                               paste0("--clump-snp-field ", opt$snp,sep=""), paste0("--clump-p1 ", opt$p1,sep=""),paste0("--clump-r2 ",opt$rsq,sep=""),paste0("--clump-p2 ",opt$p2,sep=""),paste0("--clump-kb ",opt$kb,sep=""),paste0("--out CLUMPED_",opt$out,"_",opt$p1,"_",opt$rsq,"_",opt$kb,"_chr",i,sep="")," --silent"))
    }



if (is.null(list.files(pattern="CLUMPED*.clumped"))){
  stop("ERROR: Check sumstats colnames or SNP names")}

list_clumped<-list.files(pattern=paste0("CLUMPED_",opt$out,"_",opt$p1,"_",opt$rsq,"_",opt$kb,"_chr","*.*clumped"))
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

write.table(clump_all,paste0("index_snps_",opt$out,"_",opt$p1,"_",opt$rsq,"_",opt$kb,".txt",sep=""),col.names = F,row.names=F,append=F,quote=F)

rmv<-c(paste0("CLUMPED_",opt$out,"_",opt$p1,"_",opt$rsq,"_",opt$kb,"|(_chr)|",".clumped"))
list_extract<-gsub(rmv,"",list_clumped)

for (i in list_extract){
  bin<-midf[midf$chr==i,"archivo"]
  system2("plink",args = c(paste0("--bfile ",bin,sep=""),"--allow-no-sex --extract ",
                           paste0("index_snps_",opt$out,"_",opt$p1,"_",opt$rsq,"_",opt$kb,".txt")," --make-bed ",
                           paste0("--out extracted_snps_index_",opt$out,"_",opt$p1,"_",opt$rsq,"_",opt$kb,"_chr",i,sep="")," --silent"))}

if (is.null(list.files(pattern="extracted*.bim"))){stop("ERROR: No SNPs extracted")}

system(paste0("ls extracted_snps_index_",opt$out,"_",opt$p1,"_",opt$rsq,"_",opt$kb,"_chr*.bim | sed \"s/.bim//g\" > merge_list",opt$out,opt$p1,"_",opt$rsq,"_",opt$kb,".txt",sep=""))

system2("plink",args = c("--allow-no-sex ",paste0("--merge-list merge_list",opt$out,opt$p1,"_",opt$rsq,"_",opt$kb,".txt"),
                         "--indiv-sort 0 ","--make-bed ",paste0("--out ALLclumped_",opt$out,"_",opt$p1,"_",opt$rsq,"_",opt$kb,sep="")," --silent")) 



gwas_all<-read.table(opt$summary,header=T)

clumpeados_all<-read.table(paste0("index_snps_",opt$out,"_",opt$p1,"_",opt$rsq,"_",opt$kb,".txt",sep=""))
todo_all<-merge(gwas_all,clumpeados_all,by.y="V1",by.x=opt$snp)
write.table(todo_all,paste0("Summary_clumped_",opt$out,".txt"),row.names=F,quote=F)

