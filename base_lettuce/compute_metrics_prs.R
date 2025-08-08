library(optparse)

option_list <- list(
  make_option(c("--wd"), type = "character", default = NULL, 
              help = "Working directory", metavar = "character"),
  make_option(c("--score"), type = "character", default = NULL,
              help = "Score file prefix", metavar = "character"),
  make_option(c("--p1"), type = "double", default = 1e-05,
              help = "P-val for clumping. Default 1e-05", metavar = "double"),
  make_option(c("--rsq"), type = "double", default = 0.5,
              help = "LD-threshold for clumping. Default 0.5", metavar = "double"),
  make_option(c("--kb"), type = "integer", default = 250,
              help = "Kb distance for clumping. Default 250", metavar = "integer"),
  make_option(c("--pheno"), type = "character", default = NULL,
              help = "Pheno file (txt SEPARATED BY SPACES)", metavar = "character"),
  make_option(c("--pheno_field"), type = "character", default = NULL,
              help = "Pheno field on pheno file", metavar = "character"),
  make_option(c("--covar"), type = "character", default = NULL,
              help = "Covariates for the model separated by comma", metavar = "character"),
  make_option(c("--nrep"), type = "integer", default = 1000,
              help = "Bootstrap reps. Default 1000", metavar = "integer"),
  make_option(c("--out"), type = "character", default = "out",
              help = "Output prefix for all files (analysis name)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Input validation
if (is.null(opt$wd)) stop("Working directory (--wd) is required")
if (is.null(opt$score)) stop("Score file prefix (--score) is required")
if (is.null(opt$pheno)) stop("Phenotype file (--pheno) is required")
if (is.null(opt$pheno_field)) stop("Phenotype field (--pheno_field) is required")
if (is.null(opt$covar)) stop("Covariates (--covar) are required")

wd <- opt$wd
if (!dir.exists(wd)) stop(paste("Working directory does not exist:", wd))
setwd(wd)

# Load required libraries
if (!requireNamespace("rms", quietly = TRUE)) {
  stop("rms package is required but not installed")
}
if (!requireNamespace("boot", quietly = TRUE)) {
  stop("boot package is required but not installed")
}

library(rms)
library(boot)

if (!file.exists(opt$pheno)) stop(paste("Phenotype file does not exist:", opt$pheno))

pheno <- read.table(opt$pheno, header = TRUE, stringsAsFactors = FALSE)


pheno[["fenotipo"]] <- pheno[[opt$pheno_field]]
pheno$fenotipo <- as.factor(pheno$fenotipo)

if (length(levels(pheno$fenotipo)) != 2) {
  warning("Phenotype should be binary for logistic regression analysis")
}

score_file <- paste0(opt$score, "_", opt$p1, "_", opt$rsq, "_", opt$kb, ".profile")

prs <- read.table(score_file, header = TRUE)
pheno2 <- merge(pheno, prs, by = "IID")
pheno2$st.score <- (pheno2$SCORE - mean(pheno2$SCORE)) / sd(pheno2$SCORE)

covariates <- unlist(strsplit(opt$covar, ","))
f1 <- as.formula(paste("fenotipo ~", paste(covariates, collapse = " + ")))
f2 <- as.formula(paste("fenotipo ~", paste(c(covariates, "st.score"), collapse = " + ")))


dif_rsq <- function(data, index, formula1, formula2) {
  bt <- data[index, ]
  
  mod1 <- tryCatch({
    lrm(formula1, bt, x = TRUE, y = TRUE)
  }, error = function(e) {
    warning(paste("Model 1 fitting failed:", e$message))
    return(NULL)
  })
  
  mod2 <- tryCatch({
    lrm(formula2, bt, x = TRUE, y = TRUE)
  }, error = function(e) {
    warning(paste("Model 2 fitting failed:", e$message))
    return(NULL)
  })
  
  if (is.null(mod1) || is.null(mod2)) {
    return(NA)
  }
  
  dif <- mod2$stats["R2"] - mod1$stats["R2"]
  return(dif)
}
			   
set.seed(394855)
reps <- opt$nrep

results <- boot(data = pheno2, statistic = dif_rsq, R = reps,  formula1 = f1, formula2 = f2)

outp <- mean(results$t, na.rm = TRUE)

cis <- boot.ci(results, type = "basic")

r2_iclow <- cis$basic[4]
r2_ichigh <- cis$basic[5]

n_fail <- sum(is.na(results$t))
n_success <- sum(!is.na(results$t))

outp2 <- data.frame(
  r2_inc = outp,
  ic_low = r2_iclow,
  ic_high = r2_ichigh,
  n_bootstrap = reps,
  n_success = n_success,
  n_fail = n_fail,
  PRS = paste0(opt$p1, "_", opt$rsq, "_", opt$kb)
)

output_file <- paste0(
  opt$out, "_", opt$p1, "_", opt$rsq, "_", opt$kb,
  "_validation_bootstrap_", opt$nrep, ".txt"
)

write.table(outp2, output_file, quote = FALSE, row.names = FALSE, sep = "\t")

warnings()
