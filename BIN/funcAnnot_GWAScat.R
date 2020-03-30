# LIBRARIES
library(data.table)
library(stringr)
library(parallel)
library(lme4)
library(ggsci)
library(RColorBrewer)
library(tidyverse)
library(treemap)
library(gprofiler2)
library(GOSemSim)
library(GO.db)
library(org.Hs.eg.db)
library(pheatmap)

# FUNCTIONS
## function to find variant consequences from CADD
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

## grep all eqtl from my permutation dataset
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

## function to read and store cadd files
readcadd <- function(i, ref){
    #command for reading file
    fname <- paste("INPUTS_OTHER/CADD/chr", i, "_cadd.txt", sep="")
    f <- fread(fname, h=T, showProgress=FALSE)

    f <- f[which(f$Pos %in% ref$BP & f[, "#Chrom"] == i), ]     #keep only same data as in reference data (ukb)
    print(paste("  Done with chromosome", i))
    return(f)
}

## function to plot mapping procedure
plotMapping <- function(mapping.after){
    #set how many plots to make and the margins
    layout(matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2, byrow = TRUE))
    par(mar=c(3, 3, 3, 3))

    #color palette
    colz <- pal_lancet(palette = "lanonc")(5)
    myPalette <- brewer.pal(3, "Set2")
    col.chr <-  

    #plot 1 -- barplot of the annotation sources
    pie(table(mapping.after$source_finalGenes), col=myPalette, labels=c("Coding\n(N=4)", "eQTL\n(N=31)", "Position\n(N=295)"), lwd=1.5)
    #p <- barplot(table(mapping.after$source_finalGenes), col=colz, main="Variant annotation", ylab="Number of annotations", xlab="Type of annotation", cex.names=1.25, cex.lab=1.50, cex.axis=1.25, ylim=c(0, 300), width=0.8, cex.main=2)
    #for (i in 1:nrow(p)){ text(x=p[i], y=table(mapping.after$source_finalGenes)[i]+5, label=table(mapping.after$source_finalGenes)[i], xpd=T, cex=1.25) }

    #plot 2 -- histogram of snp-gene number
    red.na <- mapping.after[which(is.na(mapping.after$geneList)), c("SNP", "geneList")]
    red <- mapping.after[which(!is.na(mapping.after$geneList)), c("SNP", "geneList")]
    #output data frame
    df <- as.data.frame(matrix(data=NA, nrow=nrow(mapping.after), ncol=2))
    colnames(df) <- c("SNP", "n_genes")
    for (i in 1:nrow(red)){
        #split geneList column
        n.gens <- strsplit(as.character(red$geneList[i]), ",")
        df[i, ] <- c(as.character(red$SNP[i]), length(n.gens[[1]]))
    }
    for (j in 1:nrow(red.na)){ df[(i+j),] <- c(as.character(red.na$SNP[j]), 0) }
    #get max
    df$n_genes <- as.numeric(df$n_genes)
    mx <- max(df$n_genes)
    #then the plot
    par(mar=c(5, 5, 2, 2))
    hist(as.numeric(df$n_genes), breaks=mx, xaxt="none", xlab="Genes per variant", ylab="Frequency", main="", cex.axis=1.25, col=colz, cex.lab=1.50, cex.main=2)
    axis(side=1, at=seq(0.5, mx, 1), labels=seq(1, mx), tick=T, cex.axis=1.25)

    #plot 3 --  histogram of distribution per chromosome
    colorz <- colorRampPalette(c("red", "orange", "dark green", "blue"))
    colorz.chr <- colorz(length(table(mapping.after$CHR)))
    barplot(table(mapping.after$CHR), col=colorz.chr, xlab="Chromosome", ylab="Genes per chromosome", cex.lab=1.50, main="", cex.axis=1.25, cex.main=2)

    return(df)
}

## main function to guide the functional annotation of variants
AnnotateMe <- function(snp_list, genes, gtex, mapping.ens){
    #run functional annotation with cadd
    chroms <- unique(snp_list$CHR)
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
    
    #finally gene-set analysis
    print("##### Performing gene-set overlap analysis")
    library(gprofiler2)
    res <- gost(geneList, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = FALSE, exclude_iea = TRUE,
        measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 1, correction_method = "fdr",
        domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = c("GO:BP"))

    return(list(res, final.annot, geneList))

}

## function to annotate using GWAS catalog
GWAScat <- function(annot, geneList){
    # read gwas catalog
    gwas <- fread("INPUTS_OTHER/GWAS_catalog_20200310.txt", h=T)

    # do 2 ways of merging to maximize results: rsid and chr:pos
    gwas$LOCUS <- paste(gwas$CHR_ID, gwas$CHR_POS, sep=":")
    annot$LOCUS <- paste(annot$CHR, annot$BP, sep=":")
    m1 <- gwas[which(gwas$SNPS %in% annot$V1), ]
    m2 <- gwas[which(gwas$LOCUS %in% annot$LOCUS), ]
    all.m <- rbind(m1, m2)
    all.m <- all.m[!duplicated(all.m), ]
    
    # separate data to be plotted
    traits <- as.data.frame(table(all.m$"MAPPED_TRAIT"))
    traits <- traits[order(-traits$Freq),]
    traits$id <- seq(1, nrow(traits))

    # for the plot, need new names
    newNames <- c("Age at Menarche", "Balding", "Bipolar disorder", "Body weight", "Chronic kidney disease",
                "Cigarettes per day", "Thyroid cancer", "Eosinophil counts", "Intelligence", "Lipid metabolism",
                "Lipoprotein metabolism", "Loneliness", "Mathematical ability", "Mean corpuscolar hemoglobin",
                "Metabolic syndrome", "Nicotine dependence", "Obesity", "Pulse pressure", "Risk-taking behaviour", 
                "Schizofrenia", "Self-reported education", "Serum carcinoembryionic antigen", "Smoking", "Thyroid peroxidase level", 
                "Total cholesterol", "Venous thromboembolism", "Vital capacity", "Waist circumference", "Waist-hip ratio", "Height", 
                "Cardiovascular disease", "Bone mineral density", "Smoking", "Coronary Artery Disease", "Systolic blood pressure", "Parental longevity")

    # run function for lolliplot
    pdf("RESULTS/PRS_ANNOTATION/gwas_cat_snps_overlap.pdf", height=6, width=6)
    lolliPlot_snp(traits, all.m)
    dev.off()

    # maybe also check genes directly
    gwas.sb <- gwas[, c("REPORTED GENE(S)", "MAPPED_GENE", "MAPPED_TRAIT")]
    gwas.sb <- as.matrix(gwas.sb)
    intervals <- ceiling(seq(nrow(gwas.sb)/10, nrow(gwas.sb), nrow(gwas.sb)/10))
    tmp_f <- function(i, gwas.sb, intervals){
        # update
        if (i %in% intervals){ print(paste("  Performed", grep(i, intervals)*10, "% of the lines.."), sep="") }
        # isolate line
        l <- gwas.sb[i, ]
        all.g <- c(unlist(strsplit(l[1], ",")), unlist(strsplit(l[2], ",")))
        all.g <- all.g[!duplicated(all.g)]
        if (TRUE %in% table(all.g %in% geneList)){ 
            tmp.mt <- data.frame(gene = all.g, trait = rep(l[3], length(all.g)))
            return(tmp.mt) 
        } else {
            tmp.mt <- data.frame(gene = NA, trait = NA)
            return(tmp.mt)
        }
    }
    res_genes <- mclapply(1:nrow(gwas.sb), tmp_f, gwas.sb=gwas.sb, intervals=intervals, mc.cores=3)     # annotate every row (~200k rows)
    all.genes <- as.data.frame(rbindlist(res_genes))       # merge results
    all.genes <- all.genes[!is.na(all.genes), ]       # exclude NAs
    overlapping.genes <- all.genes[which(all.genes$gene %in% geneList),]        # find the overlapping genes with my gene list
    overlapping.genes$gene <- as.character(overlapping.genes$gene)      # need to change these and put as character
    overlapping.genes$trait <- as.character(overlapping.genes$trait)
    overlapping.genes <- overlapping.genes[!duplicated(overlapping.genes), ]    # exclude duplicated rows

    # sampling
    n.samp <- 1000
    samp.res <- lapply(1:n.samp, sampleAnnot, overlapping.genes=overlapping.genes, annot=annot)       #sample 1000 times 1 gene per snp and do the overlap 1000 times
    options(warn=-1)
    xxx <- Reduce(function(x, y) merge(x, y, by="Var1", all.x=T, all.y=T), samp.res)
    options(warn=0)
    xxx$SUM <- rowSums(xxx[, 2:ncol(xxx)], na.rm=T)
    xxx <- xxx[, c("Var1", "SUM")]
    xxx$MEAN <- xxx$SUM/n.samp
    xxx <- xxx[order(-xxx$MEAN),]
    xxx.rec <- xxx[which(xxx$MEAN >= 1),]    
    colz <- brewer.pal(8, "Set1")
    # finally lolliplot for genes
    pdf("RESULTS/PRS_ANNOTATION/gwas_cat_genes_overlap.pdf", height=6, width=6)
    lolliPlot_gene(xxx, overlapping.genes, geneList)
    dev.off()


    l <- list(xxx, traits)
    return(xxx)
}

## function to make lolli-plot for snps
lolliPlot_snp <- function(traits, all.m){
    # will plot top 10
    traits <- traits[order(-traits$Freq), ]
    sb.tr <- head(traits, 10)
    # set new names
    nam <- c("Serum metabolite\nmeasurement", "Low density lipoprotein\ncholesterol measurement", "Blood protein\nmeasurement", "Coronary Artery\nDisease", "Total cholesterol\nmeasurement", "Physical activity\nBMI", "Parental longevity", "Systolic blood\npressure", "Smoking", "Heel bone mineral\ndensity")
    # graphical parameters
    par(mar=c(5, 11, 4, 3))
    # basic empty plot
    plot(0, 0, pch=16, col="white", xlim=c(0, 40), ylim=c(1, 10), xlab="Overlapping SNPs", ylab="", cex.axis=1.25, cex.lab=1.50, yaxt="none", bty='n')
    # add grid
    for (i in seq(0, 40, 4)){ abline(v=i, lwd=0.4, col="grey80") }
    for (i in seq(1, 10)){ abline(h=i, lwd=0.4, col="grey80") }
    # set up colors
    colz <- colorRampPalette(c("red", "orange", "dark green", "blue"))
    colorz <- colz(nrow(sb.tr))
    # add lollipops
    c <- 1
    for (i in seq(10, 1)){
        segments(x0=0, y0=i, x1=sb.tr$Freq[c], y1=i, lwd=3, col=colorz[c])
        points(x=sb.tr$Freq[c], y=i, pch=16, col=colorz[c], cex=3)
        text(x=-1, y=i, labels=nam[c], cex=1, adj=1, xpd=T)
        c <- c+1
    }   
}

## function to make lolli-plot for genes
lolliPlot_gene <- function(xxx, overlapping.genes, geneList){
    # will plot top 10
    sb.tr <- head(xxx, 10)
    # set new names
    nam <- c("Serum IgG glycosilation\nmeasurement", "Obesity", "Crohn's disease", "Celiac disease", "High density lipoprotein\ncholesterol measurement", "Total cholesterol\nmeasurement", "Bipolar disorder", "Low density lipoprotein\ncholesterol measurement", "Mental disorder", "Height")
    # graphical parameters
    par(mar=c(5, 11, 4, 3))
    # basic empty plot
    plot(0, 0, pch=16, col="white", xlim=c(0, 15), ylim=c(1, 10), xlab="Mean Overlapping Genes", ylab="", cex.axis=1.25, cex.lab=1.50, yaxt="none", bty='n', xaxt='none')
    # add grid
    for (i in seq(0, 15, 1.5)){ abline(v=i, lwd=0.4, col="grey80") }
    for (i in seq(1, 10)){ segments(x0=0, y0=i, x1=15, y1=i, lwd=0.4, col="grey80") }
    # set up colors
    colz <- colorRampPalette(c("grey80", "orange", "deepskyblue3", "darkolivegreen3"))
    colorz <- colz(nrow(sb.tr))
    # axis
    axis(side=1, at=seq(0, 15, 5), label=seq(0, 15, 5), cex.lab=1.25, cex.axis=1.25)
    # add lollipops
    c <- 1
    for (i in seq(10, 1)){
        segments(x0=0, y0=i, x1=sb.tr$MEAN[c], y1=i, lwd=3, col=colorz[c])
        points(x=sb.tr$MEAN[c], y=i, pch=16, col=colorz[c], cex=3)
        text(x=-1, y=i, labels=nam[c], cex=1, adj=1, xpd=T)
        # grep genes
        gg <- overlapping.genes[which(overlapping.genes$trait == sb.tr$Var1[c]), "gene"]
        print(gg[which(gg %in% geneList)])
        # put genes into a rectangle in the plot -- in photoshop
        #rect(xleft=15, ybottom=i-0.40, xright=22.5, ytop=i+0.40, xpd=T, col=alpha("grey80", 0.15), border=NA)
        c <- c+1
    }   
}

## function to sample for multiple genes associated with each variant
sampleAnnot <- function(i, overlapping.genes, annot){
    all.genes <- strsplit(annot$geneList, ",")
    tmp_f <- function(x, all.genes){ tmp <- all.genes[[x]]; g <- sample(x=tmp, size=1); return(g)}
    gset <- lapply(1:length(all.genes), tmp_f, all.genes=all.genes)
    gset.clean <- unlist(gset)
    overl <- overlapping.genes[which(gset.clean %in% overlapping.genes$gene),]
    df <- as.data.frame(table(overl$trait))     # make table of the trait variable
    df <- df[order(-df$Freq),]      # order the data frame according to occurrence

    return(df)
}

## function to create n sampling dsets
sampling <- function(i, mapping){
    all.genes <- strsplit(mapping$geneList, ",")
    tmp_f <- function(x, all.genes){ tmp <- all.genes[[x]]; g <- sample(x=tmp, size=1); return(g)}
    gset <- lapply(1:length(all.genes), tmp_f, all.genes=all.genes)
    gset.clean <- unlist(gset)
    return(gset.clean)
}

## function to perform overlap analysis
overlapAnalysis <- function(i, gene_list){
    g <- gene_list[[i]]

    res <- gost(g, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = FALSE, exclude_iea = TRUE,
        measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 1, correction_method = "fdr",
        domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = c("GO:BP"))
    return(res$result)
}

## function to merge sampling procedure for the functional annotation
mergeSampling <- function(enrich.res){
    all.enrich <- do.call("rbind", enrich.res)      #put all results together
    sign.only <- all.enrich[which(all.enrich$p_value <= 0.05),]         #isolate only significant findings
    sign.table <- as.data.frame(table(sign.only$term_name))
    info <- sign.only[, c("term_name", "term_id")]
    sign.table <- sign.table[order(-sign.table$Freq),]
    sign.table <- merge(sign.table, info, by.x="Var1", by.y="term_name")

    # also try with averaging pvalues
    terms <- all.enrich[, c("term_name", "term_id")]
    terms <- terms[!duplicated(terms$term_name),]
    avgP <- function(i, terms, all.enrich){ df <- data.frame(term_name = terms$term_name[i], term_id = terms$term_id[i], avgP = mean(all.enrich$p_value[which(all.enrich$term_name == terms$term_name[i])])); return(df) }
    terms.avg <- mclapply(1:nrow(terms), avgP, terms=terms, all.enrich=all.enrich, mc.cores=4)
    all.terms.avg <- rbindlist(terms.avg)
    all.terms.avg <- all.terms.avg[order(all.terms.avg$avgP),]
    all.terms.avg$log10P <- -log10(all.terms.avg$avgP)
    pdf("RESULTS/PRS_ANNOTATION/treemap_sampling_PRS_genes.pdf", height=10, width=10)
    treemap(all.terms.avg[which(all.terms.avg$avgP <= 0.05),], index="term_name", vSize="log10P", type="index")
    dev.off()

    l <- list(sign.table, all.terms.avg)
    return(l)
}

## function to do clustering of GO based on semantic similarity
semSim <- function(avg_pvalues){
    hsGO <- godata('org.Hs.eg.db', ont="BP")        #set database
    goSim("GO:0007399", "GO:0007399", semData=hsGO, measure="Wang") # example for two terms

    # extract only go terms significant
    onlyGO <- avg_pvalues[grep("GO:", avg_pvalues$term_id), ]
    onlyGO_sign <- onlyGO[which(onlyGO$avgP <= 0.05),]
    go1 = as.character(onlyGO_sign$term_id)
    go2 = as.character(onlyGO_sign$term_id)
    mt <- mgoSim(go1, go2, semData=hsGO, measure="Jiang", combine=NULL)        # example for two vectors of terms
    # the idea is to merge similar go terms into the same category -- based on semantic similarity
    # use REVIGO (website) to input the list of significant GO and their pvalue
    # get a csv as output -- save the link of the output and put it in RESULTS/PRS_ANNOTATION/
    clusters <- read.table("RESULTS/PRS_ANNOTATION/REVIGO_lin_small.csv", h=T, sep=",")
    clusters$clus_n <- NA
    clusters$col <- NA
    n = 1
    for (x in 1:nrow(clusters)){ 
        if (clusters$eliminated[x] == 0){ 
            clusters$clus_n[x] <- n
            n <- n+1
        } else {
            clusters$clus_n[x] <- n-1
        }  
    }
    annot <- data.frame(clusters = clusters$clus_n)
    rownames(annot) <- clusters$term_ID
    annot <- annot[order(annot$clusters),]
    all <- merge(mt, annot, by="row.names")
    all <- all[order(all$clusters),]
    annot <- all[, c("Row.names", "clusters")]
    all[, "clusters"] <- NULL
    rownames(all) <- all$Row.names
    all[, "Row.names"] <- NULL
    rownames(annot) <- annot$Row.names
    annot$Row.names <- NULL
    annot$clusters <- as.character(annot$clusters)
    
    colz <- colorRampPalette(grDevices::rainbow(length(table(annot$clusters))))
    mycolors <- colz(length(unique(annot$clusters)))
    names(mycolors) <- unique(annot$clusters)
    mycolors <- list(annot = mycolors)
    
    pdf("RESULTS/PRS_ANNOTATION/corrplot_significant_terms_GO_jiang.pdf", height=12, width=12)
    pheatmap(mt)
    dev.off()
}

# MAIN
## read main snp file and do the necessary adjustments
snp_list <- read.table("RESULTS/PRS/variants_5e-05.txt", h=F)     #read input list
snp_list <- as.data.frame(str_split_fixed(snp_list$V1, "_", 2))
ref <- fread("INPUTS_OTHER/GWAS_data_Timmers_cleaned.txt", h=T)
snp_list <- merge(snp_list, ref, by.x="V1", by.y="SNP")
## exclude apoe
snp_list <- snp_list[which(snp_list$BP != 45411941 & snp_list$V1 != "rs7412"),]

## other files needed
genes <- fread("INPUTS_OTHER/NCBI37.3.gene.loc", h=F, stringsAsFactors=F)
colnames(genes) <- c("gene_id", "chr", "start_tss", "stop_tss", "strand", "gene_name")
gtex <- fread("INPUTS_OTHER/Whole_Blood.v7.signif_variant_gene_pairs.txt.gz")
mapping.ens <- fread("INPUTS_OTHER/Ensemble_to_GeneName.txt")
colnames(mapping.ens) <- c("ensemble", "gene")

## manage cadd files
cadd <- mclapply(1:22, readcadd, ref=ref, mc.cores=3)

## main function for functional annotation
res <- AnnotateMe(snp_list, genes, gtex, mapping.ens)
## split the different outputs
enrich.res <- res[[1]]
enrich <- enrich.res[[1]]
geneUse <- enrich.res[[2]]
annot <- res[[2]]
geneList <- res[[3]]

## write some outputs -- geneList and variant annotation
write.table(geneList, "RESULTS/VARIANT_GENE_MAPPING/geneList_PRS.txt", quote=F, row.names=F, col.names=F)
write.table(annot, "RESULTS/VARIANT_GENE_MAPPING/PRS_variants_annotation.txt", quote=F, row.names=F, sep="\t")

## plot mapping characteristics
colnames(annot)[1] <- "SNP"
pdf("RESULTS/PRS_ANNOTATION/reported_snp_gene_mapping.pdf", height=10, width=10)
genesPerSNP <- plotMapping(kk)
dev.off()

## annotate using GWAS catalog
gwas.cat.res <- GWAScat(annot, geneList)
annot.perGene <- gwas.cat.res[[1]]
annot.perSNP <- gwas.cat.res[[2]]

## sampling to confirm functional annotation
n.sampl <- 1000
gene_sampling_dsets <- mclapply(1:n.sampl, sampling, mapping=annot, mc.cores=3)
## finally gene-set analysis for all sampling dsets
enrich.res <- mclapply(1:1000, overlapAnalysis, gene_list=gene_sampling_dsets, mc.cores=4)
sampling.res <- mergeSampling(enrich.res)
sign.table <- sampling.res[[1]]
avg_pvalues <- sampling.res[[2]]
## barplot of the top 20 terms enriched
qq <- avg_pvalues
qq <- qq[order(qq$avgP),]
qq <- avg_pvalues[1:20, ]

## save workspace
save.image("RESULTS/PRS_ANNOTATION/R_workspace_sampling_annotation.RData")

##################################################

