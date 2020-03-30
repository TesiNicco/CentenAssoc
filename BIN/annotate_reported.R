# LIBRARIES
library(ggsci)
library(data.table)
library(stringr)
library(parallel)
library(lme4)

# FUNCTIONS 
## function to read and store cadd files
readcadd <- function(i, ref){
    #command for reading file
    fname <- paste("INPUTS_OTHER/CADD/chr", i, "_cadd.txt", sep="")
    f <- fread(fname, h=T, showProgress=FALSE)

    f <- f[which(f$Pos %in% ref$BP & f[, "#Chrom"] == i), ]     #keep only same data as in reference data (ukb)
    print(paste("  Done with chromosome", i))
    return(f)
}

## alternative to the previous function to read instead of grepping at the chromsome level
SNP_conseq_CHROM_generic <- function(i, df, cadd){  
    #make subset in the input variants according to chromosome
    sb <- df[which(df$CHR == i),]

    #take subset from bigger file
    d <- cadd[[i]]

    #in case of duplicates for the same variant, isolate them and take the one with non-missing GeneName field
    sb$snp_conseq <- NA
    sb$snp_conseq_gene <- NA
    for (x in 1:nrow(sb)){
        all <- d[which(d$Pos == sb$BP[x]), ]
        info <- paste0(unique(all$Consequence[!is.na(all$Consequence)]), collapse=",")
        g <- paste0(unique(all$GeneName[!is.na(all$GeneName)]), collapse=",")
        if (info != "") { sb$snp_conseq[x] <- info } else { sb$snp_conseq[x] <- NA }
        if (g != "") { sb$snp_conseq_gene[x] <- g } else { sb$snp_conseq_gene[x] <- NA}
    }
    print(paste("  Done with chromosome", i))

    return(sb)
}

## function to put information about snp consequences with the permutation sets
mergeInfo_generic <- function(info.sb){
    #assign flag for coding snp --> culprit gene known
    info.sb$coding_snp <- NA
    for (j in 1:nrow(info.sb)){
        info <- strsplit(info.sb$snp_conseq[j], ",")[[1]]
        culprit <- unique(info %in% c("NON_SYNONYMOUS", "SYNONYMOUS", "STOP_GAINED", "MISSENSE"))    #if the variant is not missense, synonymous, stop_gained
        if (TRUE %in% culprit){
            info.sb$coding_snp[j] <- "yes"
        }
    }
    return(info.sb)
}

## grep all eqtl from my permutation dataset
GTEx_me_generic <- function(mapping, gtex, mapping.ens){
    tmp <- str_split_fixed(gtex$variant_id, "_", 5)
    gtex$newcode <- paste(tmp[, 1], "_", tmp[, 2], "_", sep="")

    mapping$eqtl_blood <- NA
    mapping$code_gtex <- paste(mapping$CHR, "_", mapping$BP, "_", sep="")
    gtex.sb <- gtex[which(gtex$newcode %in% mapping$code_gtex), ]

    for (j in 1:nrow(mapping)){
        if (is.na(mapping$coding_snp[j])){ 
            #look for eqtl
            eqtl <- gtex.sb[which(gtex.sb$newcode == mapping$code_gtex[j]), ]
            if (nrow(eqtl) >0){ #if there are eqtl, save results and stop here
                #grep gene name
                id1 <- str_split_fixed(eqtl$gene_id, "\\.", 2)
                gene_name <- mapping.ens[which(mapping.ens$ensemble %in% id1[, 1]), "gene"]
                if (nrow(gene_name) >1){ 
                    mapping$eqtl_blood[j] <- paste0(gene_name$gene, collapse=",")
                } else if (nrow(gene_name) ==1) { 
                    mapping$eqtl_blood[j] <- as.character(gene_name[, 1]) 
                } 
            } else {
                mapping$eqtl_blood[j] <- NA
            }
        }
    }
    return(mapping)
}

## function to run positional mapping for the variants with no annotation from consequences and eqtl
positionalAnn_generic <- function(mapping, genes){
    mapping$positional_mapping <- NA        #add column to be filled
    tmp.exc <- mapping[!(is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood)), ]      #identify snps to not look at
    mapping <- mapping[is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood), ]     #make clean dataset
    
    for (j in 1:nrow(mapping)){
        w = 500000   #set window threshold -- this value will be doubled (upstream and downstream)
        sb.gene <- genes[which(genes$chr == mapping$CHR[j]),]       #restrict to chromosome of use
        sb.gene <- sb.gene[which((abs(mapping$BP[j] - sb.gene$start_tss <= w)) | (abs(mapping$BP[j] - sb.gene$stop_tss <= w))), ]       #restrict within 1MB region
        sb.gene$distance_start <- abs(mapping$BP[j] - sb.gene$start_tss)        #calculate distance from start tss
        sb.gene$distance_stop <- abs(mapping$BP[j] - sb.gene$stop_tss)      #calculate distance from stop tss
        sb.gene$distance_start[which(mapping$BP[j] >= sb.gene$start_tss & mapping$BP[j] <= sb.gene$stop_tss)] <- 0      #also assign 0 for intronic positions
        sb.gene$min_distance <- apply(X=sb.gene[, c(7, 8)], MARGIN=1, FUN=min)      #function to get the minimum distance between t-start-s and t-stop-site 
        sb.gene$code <- apply(X=sb.gene[, c(9)], MARGIN=1, FUN=function(x){return(round(x/(w/10), 0))})       #assign code (intron=0, <50kb=1, <100kb=2, etc)
        
        if (nrow(sb.gene) > 0){  #if there are hits, save results and stop here
            mapping$positional_mapping[j] <- paste0(sb.gene$gene_name[which(sb.gene$code == min(sb.gene$code))], collapse=",")
        }
    }
    mapping <- rbind(mapping, tmp.exc)      #re-join dataset
    return(mapping)
}

## function to clean file after mapping procedure and output the list of genes
getGeneList_mod_generic <- function(mapping){
    #put information about which source will be used to get genes: priority is: 1-coding_snp -- 2-eqtl -- 3-positional
    mapping$source_finalGenes <- NA
    mapping$geneList <- NA

    #assign source and gene list for coding variants
    mapping$source_finalGenes[which(!is.na(mapping$coding_snp))] <- "coding"
    mapping$geneList[which(!is.na(mapping$coding_snp))] <- mapping$snp_conseq_gene[which(!is.na(mapping$coding_snp))]
    
    #assign source for eqtl hits
    mapping$source_finalGenes[which(!is.na(mapping$eqtl_blood))] <- "eqtl"
    mapping$geneList[which(!is.na(mapping$eqtl_blood))] <- mapping$eqtl_blood[which(!is.na(mapping$eqtl_blood))]

    #assign source for eqtl hits
    mapping$source_finalGenes[which(is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood))] <- "positional"
    mapping$geneList[which(is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood))] <- mapping$positional_mapping[which(is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood))]

    #assign finally those that did not map anywhere
    mapping$source_finalGenes[which(is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood) & is.na(mapping$positional_mapping))] <- "missing"

    tmp1 <- mapping$geneList[which(!is.na(mapping$geneList))]
    tmp2 <- strsplit(tmp1, ",")
    all.genes <- c()
    for (i in 1:length(tmp2)) { all.genes <- c(all.genes, tmp2[[i]]) }
    all.genes <- all.genes[!is.na(all.genes)]
    all.genes <- all.genes[which(!(all.genes %in% c("NA", "character(0)")))]
    all.genes <- all.genes[!duplicated(all.genes)]
    all.genes <- str_split_fixed(all.genes, "\\.", 2) 
    all.genes <- all.genes[, 1]

    #merge results
    l = list(all.genes, mapping)
    
    return(l)
}

## function to guide the functional annotation calling all the different functions in sequence
AnnotateMe <- function(snp_list, genes, gtex, mapping.ens){
    #run functional annotation with cadd
    chroms <- unique(snp_list$CHR)
    print("  Retrieving variant consequences")
    out <- mclapply(chroms, SNP_conseq_CHROM_generic, df=snp_list, cadd=cadd, mc.cores=3)
    #merge results together
    snps.conseq <- out[[1]]
    for (i in 2:length(out)){ snps.conseq <- rbind(snps.conseq, out[[i]]) }
    #manually change APOE4 cause for some reasons it is missing
    snps.conseq[which(snps.conseq$BP == 45411941), "snp_conseq"] <- "MISSENSE"
    snps.conseq[which(snps.conseq$BP == 45411941), "snp_conseq_gene"] <- "APOE"

    #need to re-couple with the permutation sets now
    out.annot <- mergeInfo_generic(info.sb=snps.conseq)

    #now run the GTEx and positional annotations
    print("##### Running GTEx annotation for all permutation datasets")
    out.gtex <- GTEx_me_generic(mapping=out.annot, gtex=gtex, mapping.ens=mapping.ens)

    #for the remaining, need to run the positional mapping
    print("##### Doing positional annotation for all permutation datasets")
    full.annot <- positionalAnn_generic(mapping=out.gtex, genes=genes)

    #extract gene list 
    print("##### Extracting unambiguous gene list for all permutation datasets")
    res.clean <- getGeneList_mod_generic(mapping=full.annot)
    geneList <- res.clean[[1]]
    final.annot <- res.clean[[2]]
    
    return(list(geneList, final.annot))
}

## function to plot variant-gene mapping
plotMapping <- function(mapping.after){
    #set how many plots to make
    par(mfrow=c(3, 1), mar=c(4, 5, 4, 5))

    #color palette
    colz <- pal_lancet(palette = "lanonc")(5)

    #plot 1 -- barplot of the annotation sources
    p <- barplot(table(mapping.after$source_finalGenes), col=colz, main="Variant annotation", ylab="Number of annotations", xlab="Type of annotation", cex.names=1.25, cex.lab=1.50, cex.axis=1.25, ylim=c(0, 30), width=0.8, cex.main=2)
    for (i in 1:nrow(p)){ text(x=p[i], y=table(mapping.after$source_finalGenes)[i]+5, label=table(mapping.after$source_finalGenes)[i], xpd=T, cex=1.25) }

    #plot 2 -- histogram of snp-gene number
    red.na <- mapping.after[which(is.na(mapping.after$geneList)), c("SNP", "geneList")]
    red <- mapping.after[which(!is.na(mapping.after$geneList)), c("SNP", "geneList")]
    #output data frame
    df <- as.data.frame(matrix(data=NA, nrow=nrow(mapping.after), ncol=2))
    colnames(df) <- c("SNP", "n_genes")
    for (i in 1:nrow(red)){
        #split geneList column
        n.gens <- strsplit(red$geneList[i], ",")
        df[i, ] <- c(as.character(red$SNP[i]), length(n.gens[[1]]))
    }
    for (j in 1:nrow(red.na)){ df[(i+j),] <- c(as.character(red.na$SNP[j]), 0) }
    #get max
    df$n_genes <- as.numeric(df$n_genes)
    mx <- max(df$n_genes)
    #then the plot
    hist(as.numeric(df$n_genes), breaks=mx, xaxt="none", xlab="Number of genes per variant", ylab="Frequency", main="Number of genes associated with variant", cex.axis=1.25, col=colz, cex.lab=1.50, cex.main=2)
    axis(side=1, at=seq(0.5, mx, 1), labels=seq(1, mx), tick=T, cex.axis=1.25)

    #plot 3 --  histogram of distribution per chromosome
    barplot(table(mapping.after$CHR), col=colz, xlab="Chromosome", ylab="Number of genes per chromosome", cex.lab=1.50, main="Number of genes per chromosome", cex.axis=1.25, cex.main=2)

    return(df)
}

## function to extract MAGMA coordinates
extractMagmaCoord <- function(geneList, genes){
	#find overlapping genes
	sb <- genes[which(genes$gene_name %in% geneList),]
	
	#missing
	miss <- geneList[which(!(geneList %in% genes$gene_name))]

	return(list(sb, miss))
}

#######################################
#######################################
# MAIN
## read input snps
d <- fread("RESULTS/SINGLE_VARIANT/reported_variants_info.txt", h=F)
colnames(d) <- c("CHR", "SNP", "BO", "BP", "REF", "ALT")

## read additional required files
genes <- fread("INPUTS_OTHER/NCBI37.3.gene.loc", h=F, stringsAsFactors=F)
colnames(genes) <- c("gene_id", "chr", "start_tss", "stop_tss", "strand", "gene_name")
gtex <- fread("INPUTS_OTHER/Whole_Blood.v7.signif_variant_gene_pairs.txt.gz")
mapping.ens <- fread("INPUTS_OTHER/Ensemble_to_GeneName.txt")
colnames(mapping.ens) <- c("ensemble", "gene")

## first, read CADD scores and divide them per chromosome
print("  Storing CADD information")
cadd <- mclapply(1:22, readcadd, ref=d, mc.cores=3)

## then run the actual annotation
res <- AnnotateMe(d, genes, gtex, mapping.ens)
geneList <- res[[1]]
final.annot <- res[[2]]

#write outputs
write.table(final.annot, "RESULTS/VARIANT_GENE_MAPPING/reported_var_annotation.txt", quote=F, row.names=F, sep="\t")
write.table(geneList, "RESULTS/VARIANT_GENE_MAPPING/geneList_reported_var.txt", quote=F, row.names=F, sep="\t")

## finally make plot of the mapping procedure
pdf("RESULTS/VARIANT_GENE_MAPPING/reportedVariants_mapping.pdf", height=10, width=10)
genesPerVar <- plotMapping(final.annot)
dev.off()

## find gene coordinates for MAGMA test
res <- extractMagmaCoord(geneList, genes)
coord <- res[[1]]
miss <- res[[2]]
write.table(coord, "RESULTS/GENE_BASED/gene_coordinates.txt", quote=F, row.names=F, sep="\t")
write.table(miss, "RESULTS/GENE_BASED/gene_UNmapped.txt", quote=F, row.names=F, sep="\t", col.names=F)
