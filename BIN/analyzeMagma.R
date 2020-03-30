# LIBRARY
library(data.table)
library(ggplot2)
library(ggsci)

#read magma results -- multi models
m.no <- fread("RESULTS/GENE_BASED/chrAll_multiSNP_2kb_annot.genes.out", h=T, stringsAsFactors=F)

#read genes code
g <- read.table("RESULTS/GENE_BASED/gene_coordinates.txt", h=T, stringsAsFactors=F)
dg <- g[, c("gene_id", "gene_name")]

#functio to do staff
f <- function(d, dg, inf){

  d.m <- merge(d, dg, by.x="GENE", by.y="gene_id", all.x=T)

  d.m <- d.m[order(d.m$CHR, d.m$START),]

  if (inf == "pcr"){
    d.m$P_ADJ <- p.adjust(d.m$P, method="fdr")
  } else {
    d.m$P_ADJ <- p.adjust(d.m$P_JOINT, method="fdr")
  }

  print(d.m[which(d.m$P_ADJ <= 0.20),])

  return(d.m)
}

#run for all datasets
info.mul.no <- f(m.no, dg, "m")

#also show results for the snp-wise top only
print("### Specifically SNP-wise TOP")
info.mul.no$P_ADJ_TOP <- p.adjust(info.mul.no$P_SNPWISE_TOP1, method="fdr", n=nrow(info.mul.no))
print(info.mul.no[which(info.mul.no$P_ADJ_TOP <= 0.20),])

#read magma results -- multi models
write.table(info.mul.no, "RESULTS/GENE_BASED/chrAll_multiSNP_2kb_annot.genes.out", quote=F, row.names=F)

