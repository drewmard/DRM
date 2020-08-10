# load
library('BEDMatrix')
library(data.table)
library(parallel)

# system arguments
args = commandArgs(trailingOnly=TRUE)
genotype <- args[1]
phenotype <- args[2]
phenoName <- args[3]
output <- args[4]
num_cores <- args[5]

# read in genotype data
geno <- BEDMatrix(genotype)
if (length(args) == 7) {
  startInd <- args[6]
  endInd <- args[7]
} else if (length(args) == 6) {
  startInd <- args[6]
  endInd <- length(colnames(geno))
} else if (length(args) == 5) {
  startInd <- 1
  endInd <- length(colnames(geno))
} else {
  stop("Incorrect number of input arguments for the DRM script.",call. = FALSE)
}

# read in phenotype data
pheno <- fread(phenotype,data.table = F,stringsAsFactors = F); pheno <- subset(pheno,!duplicated(pheno$IID)); pheno[pheno[,phenoName]==-9,phenoName] <- NA 

# align genotype and phenotype files
geno_names <- unlist(lapply(strsplit(rownames(geno),'_'),function(x) {return(x[2])}))
ind <- which(geno_names %in% pheno$IID)
pheno <- pheno[match(geno_names[ind],pheno$IID),]

DeviationRegressionModel <- function(i) {
  if (i %% 100 == 0) {print(i)}
  SNP <- geno[ind,i]
  PHENO <- pheno[,phenoName]
  X <- as.factor(SNP)
  Y.i <- tapply(PHENO, X, median,na.rm=T)
  Z.ij <- abs(PHENO - Y.i[X])
  res <- summary(lm(Z.ij~SNP))$coef[2,]
}

Fit_Model <- function(start=1,p=5000) {
  df.save <- mclapply(start:p,DeviationRegressionModel,mc.cores=num_cores)
  df.save <- do.call(rbind,df.save)
  df.save <- as.data.frame(df.save)
  df.save$SNP <- colnames(geno)[start:p]
  colnames(df.save)[1:4] <- c('Estimate','Std. Error','t value','Pr(>|t|)')
  return(df.save)
}

df.results <- Fit_Model(start=startInd,p=endInd)

df.results <- df.results[,c(5,1:4)]
colnames(df.results) <- c('SNP','effect_size','se','t','pval')

print('Saving...')
fwrite(df.results,output,sep='\t',quote=F,col.names = T,row.names = F,na="NA")
