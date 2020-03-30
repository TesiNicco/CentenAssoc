# LIBRARY
library(data.table)
args = commandArgs(trailingOnly=TRUE)
library(stringr)

# MAIN
## read data
d <- fread(args[1], h=T, stringsAsFactors=F)

#estimate inflation factor for the different tests
chisq.snpW <- qchisq(1 - d$P_SNPWISE_MEAN, 1)
chisq.snpT <- qchisq(1 - d$P_SNPWISE_TOP1, 1)
chisq.comb <- qchisq(1 - d$P_JOINT, 1)

lambda.snpW = round((median(chisq.snpW) / qchisq(0.5, 1)), 2)
lambda.snpT = round((median(chisq.snpT) / qchisq(0.5, 1)), 2)
lambda.comb = round((median(chisq.comb) / qchisq(0.5, 1)), 2)

print(paste("## Inflation factor for SNP-wise Mean test = ", lambda.snpW, sep=""))
print(paste("## Inflation factor for SNP-wise Top test = ", lambda.snpT, sep=""))
print("####")
print(paste("## Inflation factor for overall combined test = ", lambda.comb, sep=""))
print("####")
print("####")
