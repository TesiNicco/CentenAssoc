library(data.table)

# read input files
d <- fread("RESULTS/SINGLE_VARIANT/chrAll_reported_variants_assoc.txt", h=T, sep="\t")
ref <- fread("INPUTS_OTHER/GWAS_data_Timmers_cleaned.txt.gz", h=T)
###

# merge files together
d.my <- merge(d, ref, by.x="ID", by.y="SNP")

# initialize direction column
d.my$DIRECTION <- NA

# loop to assign directions
for (i in 1:nrow(d.my)){
	if (d.my$A1.x[i] == d.my$A1.y[i]){
		if (d.my$BETA.x[i]*d.my$BETA.y[i] >0){d.my$DIRECTION[i] <- "correct"} else { d.my$DIRECTION[i] <- "not_correct" }
	} else {
                if (d.my$BETA.x[i]*d.my$BETA.y[i] <0){d.my$DIRECTION[i] <- "correct"} else { d.my$DIRECTION[i] <- "not_correct" }
	}
}

#get stats about directions
table(d.my$DIRECTION)
print(table(d.my$DIRECTION)/nrow(d.my))

print(binom.test(x=table(d.my$DIRECTION)[1], n=nrow(d.my), p=0.5))

x <- d.my
x$P_ADJ_fdr	<- p.adjust(x$P.x, method="fdr")
print("### Correcting p-values for multiple tests")
print("### FDR correction")
print(paste("###     ", nrow(x[which(x$P_ADJ_fdr <= 0.05),]), "/", nrow(x), " variants reached significant at FDR<5%", sep=""))
print(x$ID[which(x$P_ADJ_fdr <= 0.05)])
print(paste("###     ", nrow(x[which(x$P_ADJ_fdr <= 0.10),]), "/", nrow(x), " variants reached significant at FDR<10%", sep=""))
print(x$ID[which(x$P_ADJ_fdr <= 0.10)])
print(paste("###     ", nrow(x[which(x$P_ADJ_fdr <= 0.20),]), "/", nrow(x), " variants reached significant at FDR<20%", sep=""))
print(x$ID[which(x$P_ADJ_fdr <= 0.20)])

write.table(d.my, "RESULTS/SINGLE_VARIANT/chrAll_reportedVar_assoc_withDirection.txt", quote=F, row.names=F, sep="\t")

