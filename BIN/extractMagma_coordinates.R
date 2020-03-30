# LIBRARIES
library(data.table)
library(stringr)

# MAIN
## read gene db
db <- fread("INPUTS_OTHER/NCBI37.3.gene.loc", h=F, stringsAsFactors=F)
colnames(db) <- c("gene_id", "chrom", "bp_start", "bp_end", "strand", "gene_name")

#read my genes
my <- read.table("gene_based/variant_gene_mapping.txt", h=T, stringsAsFactors=F, sep="\t")

#take genes from db that are in my genes
sb <- db[which(db$gene_name %in% my$gene),]

#check which genes did not map
missing <- my[which(!(my$gene %in% sb$gene_name)),]
write.table(missing, "gene_based/missing_mapping_MAGMA.txt", quote=F, row.names=F, sep="\t")

#write output of mapped genes
write.table(sb, "gene_based/mapped_genes_MAGMA.txt", quote=F, row.names=F, sep="\t")

