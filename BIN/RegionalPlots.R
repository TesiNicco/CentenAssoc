#######################################################################
# SCRIPT TO CREATE FILES NEEDED TO PERFORM SINGLE VARIANT ASSOCIATION #
# FOR THE REGIONAL PLOTS. THE REGIONS WILL BE APOE, ABO, CDK2NB, GBX2 #
#######################################################################

# LIBRARIES
library(data.table)
library(stringr)

# FUNCTIONS
## function to read MAGMA snps
readMagmaSNPs <- function(magma.toplot){
    # define output list
    l <- list()
    for (i in 1:nrow(magma.toplot)){
        # get variant ID
        tmp_id <- magma.toplot$GENE[i]
        tmp_chr <- magma.toplot$CHR[i]
        # grep the variant ID in the SNP FILES
        cmd <- paste("grep -w ", tmp_id, " RESULTS/GENE_BASED/ANNOTATIONS/chr", tmp_chr, ".genes.annot", sep="")
        res_grep <- system(cmd, intern=TRUE)
        # split results
        res_grep <- unlist(strsplit(res_grep, "\t"))
        snps <- res_grep[3:length(res_grep)]
        intervals <- unlist(strsplit(res_grep[2], ":"))
        # save results
        tmp_list <- list(snps, intervals)
        l[[i]] <- tmp_list
    }
    return(l)
}

## function to prepare for each region to plot a "--extract file with variant IDs"
prepareExtract <- function(snps){
    # define output list to put associations
    all.assoc <- list()

    # main loop on regions
    for (x in 1:length(snps)){
        # select snp
        tmp <- unlist(strsplit(snps[x], ":"))
        id <- paste(tmp[1], tmp[2], sep="_")
        # read corresponding pvar file
        bim <- fread(paste("GENO_DATA/PLINK2/chr", tmp[1], "_plink2_upd.pvar", sep=""), h=T)
        bim <- bim[which(abs(bim$POS - as.numeric(tmp[2])) <= 500000),]
        # write IDs somewhere
        write.table(bim$ID, "RESULTS/REGIONAL_PLOTS/tmp.variants", quote=F, row.names=F, col.names=F)
        # run association
        cmd_assoc <- paste("plink2 --pfile GENO_DATA/PLINK2/chr", tmp[1], "_plink2_upd --extract RESULTS/REGIONAL_PLOTS/tmp.variants --covar GENO_DATA/20200123_covariates_kept_Age_gwas.txt --covar-name PC1,PC2,PC3,PC4,PC5 --glm hide-covar cols=+beta,+a1freq,+a1freqcc,+machr2 --out RESULTS/REGIONAL_PLOTS/", id, "_association", sep="")
        system(cmd_assoc)
        # read association file back
        ass <- fread(paste("RESULTS/REGIONAL_PLOTS/", id, "_association.PHENO1.glm.logistic", sep=""), h=T)
        # save associations
        all.assoc[[x]] <- ass
    }
    return(all.assoc)
}

## function to prepare for each region to plot a "--extract file with variant IDs" -- this is for the MAGMA intervals
prepareExtractMagma <- function(snps){
    # define output list to put associations
    all.assoc <- list()

    # main loop on regions
    for (x in 1:length(snps)){
        # select snp
        tmp <- unlist(strsplit(snps[x], ":"))
        id <- paste(tmp[1], tmp[2], tmp[3], sep="_")
        # read corresponding pvar file
        bim <- fread(paste("GENO_DATA/PLINK2/chr", tmp[1], "_plink2_upd.pvar", sep=""), h=T)
        bim <- bim[which(abs(bim$POS - as.numeric(tmp[2])) <500000 | abs(bim$POS - as.numeric(tmp[2])) <500000),]
        # write IDs somewhere
        write.table(bim$ID, "RESULTS/REGIONAL_PLOTS/tmp.variants", quote=F, row.names=F, col.names=F)
        # run association
        cmd_assoc <- paste("plink2 --pfile GENO_DATA/PLINK2/chr", tmp[1], "_plink2_upd --extract RESULTS/REGIONAL_PLOTS/tmp.variants --covar GENO_DATA/20200123_covariates_kept_Age_gwas.txt --covar-name PC1,PC2,PC3,PC4,PC5 --glm hide-covar cols=+beta,+a1freq,+a1freqcc,+machr2 --out RESULTS/REGIONAL_PLOTS/", id, "_association", sep="")
        system(cmd_assoc)
        # read association file back
        ass <- fread(paste("RESULTS/REGIONAL_PLOTS/", id, "_association.PHENO1.glm.logistic", sep=""), h=T)
        # save associations
        all.assoc[[x]] <- ass
    }
    return(all.assoc)
}

# MAIN
## read single-variant reported association
repor <- fread("RESULTS/SINGLE_VARIANT/chrAll_reportedVar_assoc_withDirection.txt", h=T, sep="\t")
## read single-variant haplotype
haplo <- fread("RESULTS/SINGLE_VARIANT_HAPLOTYPE/chrAll_reported_and_haplotype.txt", h=T, sep="\t")
## read MAGMA
magma <- fread("RESULTS/GENE_BASED/chrAll_multiSNP_2kb_annot.genes.out", h=T)

## get the reported variants to plot
repor$P_ADJ <- p.adjust(repor$P.x, method="fdr", n=nrow(repor))
repor.toplot <- repor[which(repor$P_ADJ <= 0.10),]
repor.toplot$LOCUS <- paste(repor.toplot$"#CHROM", repor.toplot$POS, repor.toplot$ID, sep=":")

## get the haplotype variants to plot
haplo$P_ADJ <- p.adjust(haplo$P, method="fdr", n=nrow(haplo))
haplo.toplot <- haplo[which(haplo$P_ADJ <= 0.10),]
haplo.toplot <- haplo.toplot[which(!(haplo.toplot$ID %in% repor.toplot$ID)),]
haplo.toplot$LOCUS <- paste(haplo.toplot$"#CHROM", haplo.toplot$POS, haplo.toplot$ID, sep=":")

## get the magma genes to plot
magma.toplot <- magma[which(magma$P_ADJ_TOP <= 0.10),]
## for the magma genes, need also the actual snps included in the test
info_magma <- readMagmaSNPs(magma.toplot)
snps_info_magma <- c()
for (i in 1:length(info_magma)){ tmp <- info_magma[[i]]; info <- tmp[[2]]; snps_info_magma <- c(snps_info_magma, paste(info[1], info[2], info[3], sep=":")) }
all.assoc.magma <- prepareExtractMagma(snps_info_magma)

## get all associations in the selected regions
all.assoc <- prepareExtract(c(repor.toplot$LOCUS, haplo.toplot$LOCUS))

############################################################
## PREPARATION PART IS OVER -- NOW NEED TO START PLOTTING  #
############################################################

# LIBRARIES
library(data.table)
library(stringr)
library(ggplot2)

# FUNCTIONS
## extract recombination map and give coordinates of interests as output
findRecomb <- function(pos, w, rec){
  min.p <- pos-w
  max.p <- pos+w
  
  #extract interval of interest
  recomb <- rec[which(rec$position >= min.p & rec$position <= max.p),]
  
  return(recomb)
}

## function to parse gene location file and grep only region of interest
function.dynamicGene <- function(pos, w, chr, gene.db){
  #define searching space for genes
  min.x <- pos-w
  max.x <- pos+w
  
  #find genes in interval (min.x -- max.x) -- then clean up a bit
  gene.loc.res <- subset(gene.db, gene.db$chrom == paste("chr", chr, sep=""))
  gene.loc.res$in.int <- 0
  gene.loc.res$in.int[which((gene.loc.res$txStart >= min.x) & (gene.loc.res$txEnd <= max.x))] <- 1
  gene.loc.res$in.int[which((gene.loc.res$txStart <= min.x) & (gene.loc.res$txEnd >= min.x))] <- 1
  gene.loc.res$in.int[which((gene.loc.res$txStart <= max.x) & (gene.loc.res$txEnd >= max.x))] <- 1
  genes <- gene.loc.res[which(gene.loc.res$in.int == 1),]
  genes <- genes[order(-genes$txEnd),]
  genes <- genes[!duplicated(genes$"#geneName"),]
  
  #define position in the plot
  if (nrow(genes) == 2){
    genes$y <- c(-1, -2)    
  } else {
      n = ceiling(nrow(genes)/2)
      if (n != 0){
        v <- -1
        genes$y <- NA
        for (i in 1:nrow(genes)){
          genes$y[i] <- v
          v <- v-1
          if (v == -n-1){ v = -1 }
        }
      }
  }
  return(genes)
}

## function to assign dot size -- define range
function.pointSize <- function(assoc, ld.info, range, pos, id){
  ## basic size for all dots
  assoc$size <- 1.25
  
  ## depeding on pvalue, increase dot size
  for (x in range){ assoc$size[which(-log10(assoc$P) >= x & assoc$POS %in% ld.info$BP_B)] <- x }
  
  ## also point pch
  assoc$pch <- 16

  ## specific info for reported/leading variant
  assoc$pch[which(assoc$ID == id)] <- 23
  try(if (assoc$size[which(assoc$ID == id)] < 3){ assoc$size[which(assoc$ID == id)] <- 3 }, silent = T)
  assoc$pch[which(assoc$ID == id)] <- 24
  if (assoc$size[which(assoc$ID == id)] < 3){ assoc$size[which(assoc$ID == id)] <- 3 }
  
  return(assoc)
}

## function to assign dot color -- define range
function.LDcolor <- function(ld.info, snp.info){
  ld.info$color <- "deepskyblue2"
  
  ld.info$color[which(ld.info$R2 >= 0.4)] <- "yellow"
  ld.info$color[which(ld.info$R2 >= 0.6)] <- "orange"
  ld.info$color[which(ld.info$R2 >= 0.8)] <- "red"
  
  snp.info$color <- "black"
  snp.info$color[which(snp.info$POS %in% ld.info$BP_B)] <- ld.info$color[which(ld.info$BP_B %in% snp.info$POS)]
  return(snp.info)
}

## basic function to plot
function.plot <- function(pos, assoc, ld.info, w, genes, recomb, id, g.name, letters, i){
    ## graphical parameters for margins
    par(mar=c(4, 5, 6, 6))
  
    ## define max y
    y.lim <- 13
  
    ## these are standard values for minimun and maximum recombination rates genome-wide for hg19
    chromatin.lower <- 0
    chromatin.upper <- 100
  
    ## assign dot sizes (function for this)
    range = seq(1.5, 5, 0.75)
    snp.info <- function.pointSize(assoc, ld.info, range, pos, id)

    ## independently from plot type, I need the genes that are in the window to adjust y-axis -- here it is
    if (nrow(genes) > 0){ min.y <- min(genes$y) } else { min.y <- 0 }

    ## get color for LD information
    snp.info <- function.LDcolor(ld.info, snp.info)
  
    ## background plot
    plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=1.5, xaxt='none', ylab="", ylim=c(min.y*y.lim/12, y.lim), cex.axis = 1.25, 
        pch=16, col="white", cex=2, type = "p", xaxs="i", yaxt='none', xlim=c(min(snp.info$POS), max(snp.info$POS)), cex.main=2.50, bty='n')
    
    ## add grid
    #for (x in seq(0, y.lim, (y.lim-min.y*y.lim/10)/10)){abline(h=x, lwd=0.4, col="grey80")}
    #for (x in seq(min(snp.info$POS), max(snp.info$POS), (max(snp.info$POS)-min(snp.info$POS))/10)){ segments(x0 = x, y0 = 0, x1 = x, y1 = y.lim, col = "grey80", lwd=0.4) }
    
    ## add recombination rates: for this, need to normalize between 0 and max.y the recombination rate
    recomb$norm.rate <- (y.lim - 1) * ((recomb[, 2] - chromatin.lower) / (chromatin.upper - chromatin.lower))
    y.axis.recomb <- seq(0, chromatin.upper, chromatin.upper/4)
    y.axis.norm <- (y.lim - 1) * ((y.axis.recomb - min(y.axis.recomb))/(max(y.axis.recomb) - min(y.axis.recomb)))
    points(recomb$position, recomb$norm.rate, type="l", lwd=1.5, col="darkolivegreen3")
    
    ## add significance lines and corresponding legend -- for now two at 0.05 and genome-wide 5e-8
    abline(h=-log10(0.05), lty=2, col=alpha("darkgreen", 1))
    abline(h=-log10(5e-8), lty=2, col=alpha("purple", 1))

    ## then points
    snp.info$"-log10(P-value)" <- -log10(snp.info$P)
    other.points <- snp.info[which(snp.info$color == "black"),]
    points(x = other.points$POS, y = other.points$"-log10(P-value)", xaxt='none', pch=other.points$pch, col=alpha(other.points$color, 0.8), cex=other.points$size, type = "p")
    ld.points <- snp.info[which(snp.info$color != "black"),]
    points(x = ld.points$POS, y = ld.points$"-log10(P-value)", xaxt='none', pch=ld.points$pch, col=alpha(ld.points$color, 1), cex=ld.points$size, type = "p")
    ## finally the leading snps and reported snps
    sb <- ld.points[which(ld.points$pch %in% c(23, 24)),]
    points(x = sb$POS, y = sb$"-log10(P-value)", xaxt='none', pch=sb$pch, bg=alpha(sb$color, 1), cex=sb$size, type = "p", col="black")
    
    ## manage axes
    axes <- seq(min(snp.info$POS), max(snp.info$POS), (max(snp.info$POS) - min(snp.info$POS))/7)
    axes.labels <- round(axes/1000000, 2)
    axis(side = 1, at=axes, cex.axis=1.5, labels=axes.labels)
    axes.x <- seq(0, 12, 2)
    axis(side = 2, at = axes.x, labels = axes.x, cex.axis=1.5)
    
    ## axis for recombination rates on the right
    axis(side = 4, at = y.axis.norm, labels=seq(0, 100, 25), col='darkolivegreen3', col.axis = "darkolivegreen3", cex.axis=1.5, xpd=T)
    text(x = max(snp.info$POS), y = y.lim/5*4, "Recombination rate (cM/Mb)",srt = -90, col='darkolivegreen3', xpd=T, pos = 4, offset = 4, cex=1.5, font=2)
    text(x = min(snp.info$POS), y = y.lim/3*2, "-Log10(P-value)", srt = 90, xpd=T, pos = 2, offset = 4, cex=1.5, font=2)

    ## at last the genes  
    if (min.y != 0){
        ## manage gene names
        for (g in 1:nrow(genes)){
            ## main gene line -- full transcription sequence
            segments(x0 = genes$txStart[g], y0 = genes$y[g]*y.lim/12, x1 = genes$txEnd[g], y1 = genes$y[g]*y.lim/12, lwd=3)
            #need to divide exones from introns
            start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
            end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
            exons <- cbind(start, end)
            colnames(exons) <- c("start", "end")
            exons$start <- as.numeric(as.character(exons$start))
            exons$end <- as.numeric(as.character(exons$end))
            #main loop over exons
            for (j in 1:nrow(exons)){ rect(xleft=exons$start[j], ybottom=genes$y[g]*y.lim/12-(y.lim/4/4*0.15), xright=exons$end[j], ytop = genes$y[g]*y.lim/12+(y.lim/4/4*0.15), col="grey80") }
            if (genes$"#geneName"[g] %in% g.name){
                text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]*y.lim/12 + (y.lim/4/4*0.5), 
                labels=genes$"#geneName"[g], font=4, cex=1.25, col="coral")
            } else {
                text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]*y.lim/12 + (y.lim/4/4*0.5), 
                labels=genes$"#geneName"[g], font=4, cex=0.75)
            }
        }
    }
    
    ## finally legend
    legend("topleft", bty='n', legend = c("p=0.05", "p=5e-8"), lty=c(2,2), lwd=c(2,2), col=c("darkgreen", "purple"), cex=1.25, ncol = 1, pt.cex = 2)
    legend(x = min(snp.info$POS)+((max(snp.info$POS)-min(snp.info$POS))/2.5), y = 14, bty='n', legend = c("Rep. variant"), pch = c(24), col="black", pt.bg="red", cex=1.25, ncol = 1, pt.cex = 2)
        
    ## need to put size legened and ld colors legend on the right
    start <- min(snp.info$POS) + (max(snp.info$POS)-min(snp.info$POS))/5*4
    stp <- (max(snp.info$POS) - min(snp.info$POS))/22
    rect(xleft = start, ybottom = 7.75+4, xright = start+stp, ytop = 8+4.5, col = "deepskyblue2")
    text(x = start+stp/2, y = 8.15+4, labels = "0.2", font=2, cex=0.80)
    rect(xleft = start+stp, ybottom = 7.75+4, xright = start+stp*2, ytop = 8+4.5, col = "yellow")
    text(x = start+stp*2-stp/2, y = 8.15+4, labels = "0.4", font=2, cex=0.80)
    rect(xleft = start+stp*2, ybottom = 7.75+4, xright = start+stp*3, ytop = 8+4.5, col = "orange")
    text(x = start+stp*3-stp/2, y = 8.15+4, labels = "0.6", font=2, cex=0.80)
    rect(xleft = start+stp*3, ybottom = 7.75+4, xright = start+stp*4, ytop = 8+4.5, col = "red")
    text(x = start+stp*4-stp/2, y = 8.15+4, labels = "0.8", font=2, cex=0.80)
    text(x = start+stp*2, y = 8.5+4.5, labels = "Linkage (R2)", cex=1.25)
}

## basic function to plot haplotypes
function.plot.haplo <- function(pos, assoc, ld.info, w, genes, recomb, id, g.name, letters, i, rep){
    ## graphical parameters for margins
    par(mar=c(4, 5, 6, 6))
  
    ## define max y
    y.lim <- 13
  
    ## these are standard values for minimun and maximum recombination rates genome-wide for hg19
    chromatin.lower <- 0
    chromatin.upper <- 100
  
    ## assign dot sizes (function for this)
    range = seq(1.5, 5, 0.75)
    snp.info <- function.pointSize(assoc, ld.info, range, pos, id)

    ## independently from plot type, I need the genes that are in the window to adjust y-axis -- here it is
    if (nrow(genes) > 0){ min.y <- min(genes$y) } else { min.y <- 0 }

    ## get color for LD information
    snp.info <- function.LDcolor(ld.info, snp.info)
  
    ## background plot
    plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=1.5, xaxt='none', ylab="", ylim=c(min.y*y.lim/12, y.lim), cex.axis = 1.25, 
        pch=16, col="white", cex=2, type = "p", xaxs="i", yaxt='none', xlim=c(min(snp.info$POS), max(snp.info$POS)), cex.main=2.50, bty='n')
    
    ## add grid
    #for (x in seq(0, y.lim, (y.lim-min.y*y.lim/10)/10)){abline(h=x, lwd=0.4, col="grey80")}
    #for (x in seq(min(snp.info$POS), max(snp.info$POS), (max(snp.info$POS)-min(snp.info$POS))/10)){ segments(x0 = x, y0 = 0, x1 = x, y1 = y.lim, col = "grey80", lwd=0.4) }
    
    ## add recombination rates: for this, need to normalize between 0 and max.y the recombination rate
    recomb$norm.rate <- (y.lim - 1) * ((recomb[, 2] - chromatin.lower) / (chromatin.upper - chromatin.lower))
    y.axis.recomb <- seq(0, chromatin.upper, chromatin.upper/4)
    y.axis.norm <- (y.lim - 1) * ((y.axis.recomb - min(y.axis.recomb))/(max(y.axis.recomb) - min(y.axis.recomb)))
    points(recomb$position, recomb$norm.rate, type="l", lwd=1.5, col="darkolivegreen3")
    
    ## add significance lines and corresponding legend -- for now two at 0.05 and genome-wide 5e-8
    abline(h=-log10(0.05), lty=2, col=alpha("darkgreen", 1))
    abline(h=-log10(5e-8), lty=2, col=alpha("purple", 1))

    ## then points
    snp.info$"-log10(P-value)" <- -log10(snp.info$P)
    other.points <- snp.info[which(snp.info$color == "black"),]
    points(x = other.points$POS, y = other.points$"-log10(P-value)", xaxt='none', pch=other.points$pch, col=alpha(other.points$color, 0.8), cex=other.points$size, type = "p")
    ld.points <- snp.info[which(snp.info$color != "black"),]
    sb <- ld.points[which(ld.points$ID %in% c(id, rep)),]
    ld.points <- ld.points[which(!(ld.points$ID %in% c(rep, id))),]
    points(x = ld.points$POS, y = ld.points$"-log10(P-value)", xaxt='none', pch=ld.points$pch, col=alpha(ld.points$color, 1), cex=ld.points$size, type = "p")
    ## finally the leading snps and reported snps
    sb$pch[which(sb$ID == rep)] <- 24
    sb$pch[which(sb$ID == id)] <- 23
    points(x = sb$POS, y = sb$"-log10(P-value)", xaxt='none', pch=sb$pch, bg=alpha(sb$color, 1), cex=sb$size, type = "p", col="black")
    
    ## manage axes
    axes <- seq(min(snp.info$POS), max(snp.info$POS), (max(snp.info$POS) - min(snp.info$POS))/7)
    axes.labels <- round(axes/1000000, 2)
    axis(side = 1, at=axes, cex.axis=1.5, labels=axes.labels)
    axes.x <- seq(0, 12, 2)
    axis(side = 2, at = axes.x, labels = axes.x, cex.axis=1.5)
    
    ## axis for recombination rates on the right
    axis(side = 4, at = y.axis.norm, labels=seq(0, 100, 25), col='darkolivegreen3', col.axis = "darkolivegreen3", cex.axis=1.5, xpd=T)
    text(x = max(snp.info$POS), y = y.lim/5*4, "Recombination rate (cM/Mb)",srt = -90, col='darkolivegreen3', xpd=T, pos = 4, offset = 4, cex=1.5, font=2)
    text(x = min(snp.info$POS), y = y.lim/3*2, "-Log10(P-value)", srt = 90, xpd=T, pos = 2, offset = 4, cex=1.5, font=2)

    ## at last the genes  
    if (min.y != 0){
        ## manage gene names
        for (g in 1:nrow(genes)){
            ## main gene line -- full transcription sequence
            segments(x0 = genes$txStart[g], y0 = genes$y[g]*y.lim/12, x1 = genes$txEnd[g], y1 = genes$y[g]*y.lim/12, lwd=3)
            #need to divide exones from introns
            start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
            end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
            exons <- cbind(start, end)
            colnames(exons) <- c("start", "end")
            exons$start <- as.numeric(as.character(exons$start))
            exons$end <- as.numeric(as.character(exons$end))
            #main loop over exons
            for (j in 1:nrow(exons)){ rect(xleft=exons$start[j], ybottom=genes$y[g]*y.lim/12-(y.lim/4/4*0.15), xright=exons$end[j], ytop = genes$y[g]*y.lim/12+(y.lim/4/4*0.15), col="grey80") }
            if (genes$"#geneName"[g] %in% g.name){
                text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]*y.lim/12 + (y.lim/4/4*0.5), 
                labels=genes$"#geneName"[g], font=4, cex=1.25, col="coral")
            } else {
                text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]*y.lim/12 + (y.lim/4/4*0.5), 
                labels=genes$"#geneName"[g], font=4, cex=0.75)
            }
        }
    }
    
    ## finally legend
    legend("topleft", bty='n', legend = c("p=0.05", "p=5e-8"), lty=c(2,2), lwd=c(2,2), col=c("darkgreen", "purple"), cex=1.25, ncol = 1, pt.cex = 2)
    legend(x = min(snp.info$POS)+((max(snp.info$POS)-min(snp.info$POS))/2.5), y = 14, bty='n', legend = c("Rep. variant", "Sign. haplotype"), pch = c(24, 23), col="black", pt.bg="black", cex=1.25, ncol = 1, pt.cex = 2)
        
    ## need to put size legened and ld colors legend on the right
    start <- min(snp.info$POS) + (max(snp.info$POS)-min(snp.info$POS))/5*4
    stp <- (max(snp.info$POS) - min(snp.info$POS))/22
    rect(xleft = start, ybottom = 7.75+4, xright = start+stp, ytop = 8+4.5, col = "deepskyblue2")
    text(x = start+stp/2, y = 8.15+4, labels = "0.2", font=2, cex=0.80)
    rect(xleft = start+stp, ybottom = 7.75+4, xright = start+stp*2, ytop = 8+4.5, col = "yellow")
    text(x = start+stp*2-stp/2, y = 8.15+4, labels = "0.4", font=2, cex=0.80)
    rect(xleft = start+stp*2, ybottom = 7.75+4, xright = start+stp*3, ytop = 8+4.5, col = "orange")
    text(x = start+stp*3-stp/2, y = 8.15+4, labels = "0.6", font=2, cex=0.80)
    rect(xleft = start+stp*3, ybottom = 7.75+4, xright = start+stp*4, ytop = 8+4.5, col = "red")
    text(x = start+stp*4-stp/2, y = 8.15+4, labels = "0.8", font=2, cex=0.80)
    text(x = start+stp*2, y = 8.5+4.5, labels = "Linkage (R2)", cex=1.25)
}

## basic function to plot magma
function.plot.magma <- function(pos, assoc, snps, w, genes, recomb, id, g.name, letters, i, rep){
    ## graphical parameters for margins
    par(mar=c(4, 5, 6, 6))
  
    ## define max y
    y.lim <- 13
  
    ## these are standard values for minimun and maximum recombination rates genome-wide for hg19
    chromatin.lower <- 0
    chromatin.upper <- 100
  
    ## assign dot sizes (function for this)
    range = seq(1.5, 5, 0.75)
    assoc$size <- 1.25
    ## depeding on pvalue, increase dot size
    for (x in range){ assoc$size[which(-log10(assoc$P) >= x)] <- x }
    ## also point pch
    assoc$pch <- 16
    ## specific info for reported/leading variant
    assoc$pch[which(assoc$ID %in% snps)] <- 23

    ## independently from plot type, I need the genes that are in the window to adjust y-axis -- here it is
    if (nrow(genes) > 0){ min.y <- min(genes$y) } else { min.y <- 0 }

    ## get color for LD information
    assoc$color <- "black"
    assoc$color[which(assoc$pch == 23)] <- "coral"

    ## background plot
    plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=1.5, xaxt='none', ylab="", ylim=c(min.y*y.lim/12, y.lim), cex.axis = 1.25, 
        pch=16, col="white", cex=2, type = "p", xaxs="i", yaxt='none', xlim=c(min(assoc$POS), max(assoc$POS)), cex.main=2.50, bty='n')
    
    ## add recombination rates: for this, need to normalize between 0 and max.y the recombination rate
    recomb$norm.rate <- (y.lim - 1) * ((recomb[, 2] - chromatin.lower) / (chromatin.upper - chromatin.lower))
    y.axis.recomb <- seq(0, chromatin.upper, chromatin.upper/4)
    y.axis.norm <- (y.lim - 1) * ((y.axis.recomb - min(y.axis.recomb))/(max(y.axis.recomb) - min(y.axis.recomb)))
    points(recomb$position, recomb$norm.rate, type="l", lwd=1.5, col="darkolivegreen3")
    
    ## add significance lines and corresponding legend -- for now two at 0.05 and genome-wide 5e-8
    abline(h=-log10(0.05), lty=2, col=alpha("darkgreen", 1))
    abline(h=-log10(5e-8), lty=2, col=alpha("purple", 1))

    ## then points
    assoc$"-log10(P-value)" <- -log10(assoc$P)
    other.points <- assoc[which(assoc$color == "black"),]
    points(x = other.points$POS, y = other.points$"-log10(P-value)", xaxt='none', pch=other.points$pch, col=alpha(other.points$color, 0.8), cex=other.points$size, type = "p")
    ld.points <- assoc[which(assoc$color != "black"),]
    points(x = ld.points$POS, y = ld.points$"-log10(P-value)", xaxt='none', pch=ld.points$pch, col="black", bg=alpha(ld.points$color, 1), cex=ld.points$size, type = "p")
    ## finally the leading snps and reported snps
    
    ## manage axes
    axes <- seq(min(assoc$POS), max(assoc$POS), (max(assoc$POS) - min(assoc$POS))/7)
    axes.labels <- round(axes/1000000, 2)
    axis(side = 1, at=axes, cex.axis=1.5, labels=axes.labels)
    axes.x <- seq(0, 12, 2)
    axis(side = 2, at = axes.x, labels = axes.x, cex.axis=1.5)
    
    ## axis for recombination rates on the right
    axis(side = 4, at = y.axis.norm, labels=seq(0, 100, 25), col='darkolivegreen3', col.axis = "darkolivegreen3", cex.axis=1.5, xpd=T)
    text(x = max(assoc$POS), y = y.lim/5*4, "Recombination rate (cM/Mb)",srt = -90, col='darkolivegreen3', xpd=T, pos = 4, offset = 4, cex=1.5, font=2)
    text(x = min(assoc$POS), y = y.lim/3*2, "-Log10(P-value)", srt = 90, xpd=T, pos = 2, offset = 4, cex=1.5, font=2)

    ## at last the genes  
    if (min.y != 0){
        ## manage gene names
        for (g in 1:nrow(genes)){
            ## main gene line -- full transcription sequence
            segments(x0 = genes$txStart[g], y0 = genes$y[g]*y.lim/12, x1 = genes$txEnd[g], y1 = genes$y[g]*y.lim/12, lwd=3)
            #need to divide exones from introns
            start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
            end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
            exons <- cbind(start, end)
            colnames(exons) <- c("start", "end")
            exons$start <- as.numeric(as.character(exons$start))
            exons$end <- as.numeric(as.character(exons$end))
            #main loop over exons
            for (j in 1:nrow(exons)){ rect(xleft=exons$start[j], ybottom=genes$y[g]*y.lim/12-(y.lim/4/4*0.15), xright=exons$end[j], ytop = genes$y[g]*y.lim/12+(y.lim/4/4*0.15), col="grey80") }
            if (genes$"#geneName"[g] %in% g.name){
                text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]*y.lim/12 + (y.lim/4/4*0.5), 
                labels=genes$"#geneName"[g], font=4, cex=1.25, col="coral")
            } else {
                text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]*y.lim/12 + (y.lim/4/4*0.5), 
                labels=genes$"#geneName"[g], font=4, cex=0.75)
            }
        }
    }
    
    ## finally legend
    legend("topleft", bty='n', legend = c("p=0.05", "p=5e-8"), lty=c(2,2), lwd=c(2,2), col=c("darkgreen", "purple"), cex=1.25, ncol = 1, pt.cex = 2)
    legend(x = min(assoc$POS)+((max(assoc$POS)-min(assoc$POS))/2.5), y = 14, bty='n', legend = c("Tested variants"), pch = c(23), col="black", pt.bg="coral", cex=1.25, ncol = 1, pt.cex = 2)
}

# MAIN
## set number of plots to make
n_plots <- nrow(repor.toplot) + nrow(haplo.toplot) + nrow(magma.toplot)
## set order to plots
order <- c(rep("reported", nrow(repor.toplot)), rep("haplotype", nrow(haplo.toplot)), rep("magma", nrow(magma.toplot)))
letters <- c("A", "B", "C", "D", "E", "F")
## define global windo upstream and downstram
w <- 75000

## main loop
for (i in 1:6){
    ## check what type of plot it has to be
    type <- order[i]

    ## then address in case of reported variants
    if (type == "reported"){
        ## extract variant id, position and chromosome
        id <- repor.toplot$ID[i]
        pos <- repor.toplot$POS[i]
        chr <- repor.toplot$"#CHROM"[i]

        ## calculate LD here
        cmd <- paste("plink --bfile GENO_DATA/chr", chr, "_dose --r2 yes-really --ld-window-kb 1000 --ld-snp ", id, " --out RESULTS/REGIONAL_PLOTS/", id, "_LD_info", sep="")
        system(cmd)

        ## read file back 
        ld.info <- fread(paste("RESULTS/REGIONAL_PLOTS/", id, "_LD_info.ld", sep=""), h=T)
        
        ## read genes
        gene.db <- fread("INPUTS_OTHER/hg19_geneListandPos", h=T)
        ## then restrict to region of interest
        genes <- function.dynamicGene(pos, w, chr, gene.db)

        ## read recomb rates
        rec <- fread(paste("INPUTS_OTHER/RECOMB_RATES/genetic_map_chr", chr, "_combined_b37.txt", sep=""), h=T)
        ## then restrict to region of interest
        recomb <- findRecomb(pos, w, rec)

        ## read association file
        assoc <- all.assoc[[i]]
        ## the restrict to region of interest
        assoc <- assoc[which(abs(assoc$POS - pos) <= w),]
  
        ## time to plot now
        g.name <- c("ABO", "APOE", "GBX2", "CDKN2B")
        pdf(paste("RESULTS/REGIONAL_PLOTS/PLOT_", id, "_plot.pdf", sep=""), height=8, width=9)
        function.plot(pos, assoc, ld.info, w, genes, recomb, id, g.name, letters, i)        
        dev.off()
    } else if (type == "haplotype"){
        ## extract variant id, position and chromosome
        id <- haplo.toplot$ID[i-3]
        pos <- haplo.toplot$POS[i-3]
        chr <- haplo.toplot$"#CHROM"[i-3]
        rep <- haplo.toplot$REPORTED_VAR[i-3]

        ## calculate LD here
        cmd <- paste("plink --bfile GENO_DATA/chr", chr, "_dose --r2 yes-really --ld-window-kb 1000 --ld-snp ", id, " --out RESULTS/REGIONAL_PLOTS/", id, "_LD_info", sep="")
        system(cmd)

        ## read file back 
        ld.info <- fread(paste("RESULTS/REGIONAL_PLOTS/", id, "_LD_info.ld", sep=""), h=T)
        
        ## read genes
        gene.db <- fread("INPUTS_OTHER/hg19_geneListandPos", h=T)
        ## then restrict to region of interest
        genes <- function.dynamicGene(pos, w, chr, gene.db)

        ## read recomb rates
        rec <- fread(paste("INPUTS_OTHER/RECOMB_RATES/genetic_map_chr", chr, "_combined_b37.txt", sep=""), h=T)
        ## then restrict to region of interest
        recomb <- findRecomb(pos, w, rec)

        ## read association file
        assoc <- all.assoc[[i]]
        ## the restrict to region of interest
        assoc <- assoc[which(abs(assoc$POS - pos) <= w),]
  
        ## time to plot now
        g.name <- c("ABO", "APOE", "GBX2", "CDKN2B")
        pdf("RESULTS/REGIONAL_PLOTS/GBX2_plot.pdf", height=8, width=9)
        function.plot.haplo(pos, assoc, ld.info, w, genes, recomb, id, g.name, letters, i, rep)        
        dev.off()
    } else {
        ## extract variant id, position and chromosome
        pos.start <- magma.toplot$START[i-4]
        pos.stop <- magma.toplot$STOP[i-4]
        chr <- magma.toplot$CHR[i-4]
        snps_info <- info_magma[[i-4]]
        snps <- snps_info[[1]]
        mid.pos <- ceiling(pos.start + (pos.stop-pos.start)/2)

        ## read genes
        gene.db <- fread("INPUTS_OTHER/hg19_geneListandPos", h=T)
        ## then restrict to region of interest
        genes <- function.dynamicGene(mid.pos, w, chr, gene.db)

        ## read recomb rates
        rec <- fread(paste("INPUTS_OTHER/RECOMB_RATES/genetic_map_chr", chr, "_combined_b37.txt", sep=""), h=T)
        ## then restrict to region of interest
        recomb <- findRecomb(mid.pos, w, rec)

        ## read association file
        assoc <- all.assoc.magma[[i-4]]
        ## the restrict to region of interest
        assoc <- assoc[which(abs(assoc$POS - mid.pos) <= w),]
  
        ## time to plot now
        g.name <- c("ABO", "APOE", "GBX2", "CDKN2B")
        pdf(paste("RESULTS/REGIONAL_PLOTS/", pos.start, "_plot.pdf", sep=""), height=8, width=9)
        function.plot.magma(pos, assoc, snps, w, genes, recomb, id, g.name, letters, i, rep)        
        dev.off()

    }
}

