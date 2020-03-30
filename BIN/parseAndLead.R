# LIBRARIES
library(gridExtra)
library(data.table)
library(RColorBrewer)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

# FUNCTIONS
## define a function to do the operations at single variant level, to iterate
LookAtVar <- function(i, fileList){
	#update user
	print(paste("## Working on file:", i))

	#define file of interest and open it
	fname <- fileList[i]
	d <- fread(fname, h=T, stringsAsFactors=F)

	#take all variants in LD with reported variant and write them
	all.ld.snps <- unique(d$SNP_B)
	write.table(all.ld.snps, "RESULTS/SINGLE_VARIANT_HAPLOTYPE/tmp_snps.txt", quote=F, row.names=F, col.names=F, sep="\t")

	#need to run association for these variants
	chr <- d$CHR_A[1]
	snp_a <- d$SNP_A[1]
	cmd = paste("plink2 --bfile GENO_DATA/chr", chr, "_dose --covar GENO_DATA/20200123_covariates_kept_Age_gwas.txt --covar-name PC1,PC2,PC3,PC4,PC5 --extract RESULTS/SINGLE_VARIANT_HAPLOTYPE/tmp_snps.txt --glm hide-covar cols=+beta,+a1freq,+a1freqcc,+machr2 --out RESULTS/SINGLE_VARIANT_HAPLOTYPE/", snp_a, "_LDhaplo_assoc", sep="")
	system(cmd, ignore.stdout=T)

	#read association file back
	ass <- fread(paste("RESULTS/SINGLE_VARIANT_HAPLOTYPE/", snp_a, "_LDhaplo_assoc.PHENO1.glm.logistic", sep=""), h=T)
	ass <- ass[order(ass$P), ]
	lead <- ass[1, ]

	#add few stats like reported variant and r2 value
	lead$REPORTED_VAR <- NA
	lead$R2 <- NA
        lead$INFO <- NA
	lead$REPORTED_VAR <- snp_a
	lead$R2 <- d$R2[which(d$SNP_B == lead$ID)]
	if (lead$ID == snp_a){ lead$INFO = "Same_as_reported" } else { lead$INFO = "Different_from_reported" } 

	#finally, clean
	cmd = "rm RESULTS/SINGLE_VARIANT_HAPLOTYPE/*assoc*"
	cmd2 = "rm RESULTS/SINGLE_VARIANT_HAPLOTYPE/tmp_snps.txt"
	system(cmd)
	system(cmd2)
	
	return(lead)	
}

# MAIN
## identify all files
fileList_cmd <- "ls RESULTS/SINGLE_VARIANT_HAPLOTYPE/*ld"
fileList <- system(fileList_cmd, inter=T)

#use lapply to run over all variants
all.res <- lapply(1:length(fileList), LookAtVar, fileList=fileList)
res <- all.res[[1]]
for (i in 2:length(all.res)){ res <- rbind(res, all.res[[i]]) }
table(res$INFO)

#write output table
write.table(res, "RESULTS/SINGLE_VARIANT_HAPLOTYPE/chrAll_reported_and_haplotype.txt", quote=F, row.names=F, sep="\t")
