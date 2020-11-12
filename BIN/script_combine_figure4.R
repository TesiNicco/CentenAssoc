# revigo plot
revigoPlot <- function(one.data){
  # take only terms that were not merged
  one.data <- one.data[which(one.data$eliminated == 0),]
  
  # do some adjustements
  one.data$plot_X <- as.numeric(one.data$plot_X)
  one.data$plot_Y <- as.numeric(one.data$plot_Y)
  
  ## base plot
  par(mar=c(5, 5, 2, 8))
  plot(0, 0, pch=16, col="white", xlim=c(min(one.data$plot_X)-1, max(one.data$plot_X)+1), ylim=c(min(one.data$plot_Y)-1, max(one.data$plot_Y)+1), xlab="Semantic space X", ylab="Semantic space Y", cex.lab=1.50)
  
  ## grid
  for (i in seq(min(one.data$plot_X)-1, max(one.data$plot_X)+1, (max(one.data$plot_X)-min(one.data$plot_X)+2)/10)){abline(v=i, lwd=0.4, col="grey80")}
  for (i in seq(min(one.data$plot_Y)-1, max(one.data$plot_Y)+1, (max(one.data$plot_Y)-min(one.data$plot_Y)+2)/10)){abline(h=i, lwd=0.4, col="grey80")}
  
  ## colors
  colz <- colorRampPalette(c("red", "orange", "yellow"))
  colorz <- colz(nrow(one.data))
  one.data <- one.data[order(one.data$'log10.p.value'),]
  one.data$col <- colorz
  
  ## points
  points(x = one.data$plot_X, y = one.data$plot_Y, pch=16, cex=one.data$plot_size+5, col=colorz, lwd=2)
  
  # need to add some \n in the name description
  one.data$upd_name <- NA
  for (i in 1:nrow(one.data)){
    # count the characters
    tmp.string <- as.character(one.data$description[i])
    
    # calculate in how many lines to distribute the name -- max is 23 characters (including spaces) per line
    n.lines <- ceiling(nchar(tmp.string)/25)
    
    # divide by space depending on n.lines
    tmp <- unlist(strsplit(as.character(one.data$description[i]), " "))
    if (n.lines >1){ 
      x <- split(tmp, ceiling(seq_along(tmp)/n.lines))
      for (j in 1:length(x)){ x[[j]] <- paste(x[[j]], collapse=" ") }
      # add name to data frame
      one.data$upd_name[i] <- paste(x, collapse="\n")
    } else {
      one.data$upd_name[i] <- paste(tmp, collapse=" ")
    }
  }
  
  # uppercase first letter of sentence
  one.data$upd_name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", one.data$upd_name, perl = TRUE)
  
  ## set how many data to annotate -- top n
  if (nrow(one.data) >=25){ n <- 15 } else { n <- nrow(one.data) }
  
  ## set position for annotation
  sb <- one.data[1:n, ]
  addTextLabels(sb$plot_X, sb$plot_Y, sb$upd_name, col.label="black")
  
  ## legend
  mx.x <- max(one.data$plot_X)+2
  mx.y <- max(one.data$plot_Y)+1
  step <- 0.10
  text(x = mx.x+(step), y = mx.y-(step*3), labels = "-Log10(P)", cex=1.20, xpd=T, font=1)
  gradient.rect(xleft = mx.x-(step*2), ybottom = mx.y-(step*13), xright = mx.x+(step*2), ytop = mx.y-(step*6), col=rev(colorz), nslices = length(colorz), gradient = "horizontal")
  text(x = mx.x+(step*4), y = mx.y-(step*6), labels = round(max(abs(one.data$'log10.p.value')), 1), adj = 0.5, xpd=T, font=1, cex=0.80)
  text(x = mx.x+(step*4), y = mx.y-(step*9.5), labels = round(min(abs(one.data$'log10.p.value')), 1), adj = 0.5, xpd=T, font=1, cex=0.80)
  text(x = mx.x+(step*4), y = mx.y-(step*13), labels = round(median(abs(one.data$'log10.p.value')), 1), adj = 0.5, xpd=T, font=1, cex=0.80)
  
  text(x = mx.x, y = mx.y+(step*5), cex=2, labels = "A", font=2, xpd=T)
  # ## legend for circle size
  # text(x = mx.x, y = mx.y-(step*15), labels = "Term size", cex=1.20, font=2, xpd=T)
  # points(x = mx.x, y = mx.y-(step*18), pch=16, col="grey80", cex=min(one.data$plot_size+5), xpd=T)
  # points(x = mx.x, y = mx.y-(step*23), pch=16, col="grey80", cex=median(one.data$plot_size+5), xpd=T)
  # points(x = mx.x, y = mx.y-(step*29), pch=16, col="grey80", cex=max(one.data$plot_size+5), xpd=T)
}

## function to make lolli-plot for genes
lolliPlot_gene <- function(tb, overlapping.genes, geneList){
  # will plot top 10
  sb.tr <- head(tb, 10)
  
  # initialize variable
  sb.tr$upd_name <- NA
  sb.tr$study <- NA
  
  # need to add some \n in the name description
  for (i in 1:nrow(sb.tr)){
    # count the characters
    tmp.string <- as.character(sb.tr$Var1[i])
    tmp.study <- unique(overlapping.genes[which(overlapping.genes$trait == tmp.string), "study"])
    tmp.study <- as.data.frame(str_split_fixed(tmp.study$study, " ", 2))
    sb.tr$study[i] <- nrow(tmp.study)
    
    # calculate in how many lines to distribute the name -- max is 23 characters (including spaces) per line
    n.lines <- ceiling(nchar(tmp.string)/25)
    
    # divide by space depending on n.lines
    tmp <- unlist(strsplit(as.character(sb.tr$Var1[i]), " "))
    if (n.lines >1){ 
      x <- split(tmp, ceiling(seq_along(tmp)/n.lines))
      for (j in 1:length(x)){ x[[j]] <- paste(x[[j]], collapse=" ") }
      # add name to data frame
      sb.tr$upd_name[i] <- paste(x, collapse="\n")
    } else {
      sb.tr$upd_name[i] <- paste(tmp, collapse=" ")
    }
  }
  
  # uppercase first letter of sentence
  sb.tr$upd_name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", sb.tr$upd_name, perl = TRUE)
  
  # graphical parameters
  par(mar=c(5, 8, 3, 2))
  
  # basic empty plot
  plot(0, 0, pch=16, col="white", xlim=c(0, 70), ylim=c(1, 10), xlab="Mean Overlapping Genes", ylab="", cex.axis=1.25, cex.lab=1.50, yaxt="none", bty='n')
  
  # axis
  #axis(side=1, at=seq(0, 70, ceiling(70/8)), label=seq(0, 70, ceiling(70/8)), cex.lab=1.25, cex.axis=1.25)    
  
  # add grid
  for (i in seq(0, 70, 70/8)){ abline(v=i, lwd=0.4, col="grey80") }
  for (i in seq(1, 10)){ segments(x0=0, y0=i, x1=70, y1=i, lwd=0.4, col="grey80") }
  
  # set up colors
  colz <- colorRampPalette(c("grey80", "orange", "deepskyblue3", "darkolivegreen3"))
  colorz <- colz(nrow(sb.tr))
  
  # add lollipops
  c <- 1
  for (i in seq(10, 1)){
    segments(x0=0, y0=i, x1=sb.tr$MEAN[c], y1=i, lwd=2, col=colorz[c])
    points(x=sb.tr$MEAN[c], y=i, pch=16, col=colorz[c], cex=2, xpd=T)
    text(x=-ceiling(max(sb.tr$MEAN))*0.05, y=i, labels=sb.tr$upd_name[c], cex=0.75, adj=1, xpd=T)
    c <- c+1
  }
  
  text(x = 70, y = 11.5, cex=2, labels = "C", font=2, xpd=T)
  
}

## function to make lolli-plot for snps
lolliPlot_snp <- function(traits, all.m){
  # will plot top 10 -- isolate them here
  traits <- traits[order(-traits$Freq), ]
  sb.tr <- head(traits, 10)
  
  # initialize variable
  sb.tr$upd_name <- NA
  sb.tr$study <- NA
  
  # need to add some \n in the name description
  for (i in 1:nrow(sb.tr)){
    # count the characters
    tmp.string <- as.character(sb.tr$Var1[i])
    tmp.study <- unique(all.m[which(all.m$MAPPED_TRAIT == tmp.string), c("FIRST AUTHOR")])
    tmp.study <- as.data.frame(str_split_fixed(tmp.study$'FIRST AUTHOR', " ", 2))
    sb.tr$study[i] <- nrow(tmp.study)
    
    # calculate in how many lines to distribute the name -- max is 23 characters (including spaces) per line
    n.lines <- ceiling(nchar(tmp.string)/25)
    
    # divide by space depending on n.lines
    tmp <- unlist(strsplit(as.character(sb.tr$Var1[i]), " "))
    if (n.lines >1){ 
      x <- split(tmp, ceiling(seq_along(tmp)/n.lines))
      for (j in 1:length(x)){ x[[j]] <- paste (x[[j]], collapse=" ") }
      # add name to data frame
      sb.tr$upd_name[i] <- paste(x, collapse="\n")
    } else {
      sb.tr$upd_name[i] <- paste(tmp, collapse=" ")
    }
  }
  
  # uppercase first letter of sentence
  sb.tr$upd_name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", sb.tr$upd_name, perl = TRUE)
  sb.tr[which(sb.tr$Var1 == "low density lipoprotein cholesterol measurement"), "upd_name"] <- "LDL Cholesterol\nMeasurement"
  sb.tr[which(sb.tr$Var1 == "high density lipoprotein cholesterol measurement"), "upd_name"] <- "HDL Cholesterol\nMeasurement"

  # graphical parameters
  par(mar=c(5, 8, 3, 2))
  # basic empty plot
  plot(0, 0, pch=16, col="white", xlim=c(0, 14), ylim=c(1, 10), xlab="Overlapping SNPs", ylab="", cex.axis=1.25, cex.lab=1.50, yaxt="none", bty='n')
  # add grid
  for (i in seq(0, 14, 1.4)){ abline(v=i, lwd=0.4, col="grey80") }
  for (i in seq(1, 10)){ segments(x0 = 0, y0 = i, x1 = 14, y1 = i, lwd=0.4, col="grey80") }
  # set up colors
  colz <- colorRampPalette(c("red", "orange", "dark green", "blue"))
  colorz <- colz(nrow(sb.tr))
  # add lollipops
  c <- 1
  for (i in seq(10, 1)){
    segments(x0=0, y0=i, x1=sb.tr$Freq[c], y1=i, lwd=2, col=colorz[c])
    points(x=sb.tr$Freq[c], y=i, pch=16, col=colorz[c], cex=2, xpd=T)
    text(x=-ceiling(max(sb.tr$Freq))*0.05, y=i, labels=sb.tr$upd_name[c], cex=0.75, adj=1, xpd=T)
    c <- c+1
  }
  text(x = 14, y = 11.5, cex=2, labels = "B", font=2, xpd=T)
  
}

# layout for plot
library(basicPlotteR)
library(plotrix)
load("/Users/nicco/Desktop/2k19_work/LongevityStuff/2020_longevity_paper/RESULTS/PRS_ANNOTATION/R_workspace_sampling_annotation.RData")

d <- read.table("RESULTS/RESULTS_annotateMe_PRS/revigo_out.csv", sep=",", h=T)
all.m <- read.table("RESULTS/RESULTS_annotateMe_PRS/gwas_cat_whole.txt", h=T, sep="\t")
traits <- read.table("RESULTS/RESULTS_annotateMe_PRS/snp_GWAS_cat.txt", h=T, sep="\t")
traits_genes <- read.table("RESULTS/RESULTS_annotateMe_PRS/genes_GWAS_cat.txt", h=T, sep="\t")
geneList <- read.table("RESULTS/RESULTS_annotateMe_PRS/snp_annotation_geneList.txt", h=F)
overlapping.genes <- read.table("RESULTS/RESULTS_annotateMe_PRS/overlapping_genes_GWAScat.txt", h=T, sep="\t")

#bitmap("/Users/nicco/Desktop/myPapers/ROAD_TO_THIRD_PAPER/The journals of gerontology/Rebuttal/final_files/figure_3.tiff", width = 5, height = 5, units = "in", res = 300, type ="tifflzw")
#postscript("/Users/nicco/Desktop/myPapers/ROAD_TO_THIRD_PAPER/The journals of gerontology/Rebuttal/final_files/figure_3.eps", width = 5, height = 5, horizontal = F, paper = "special", onefile = T)
#png("/Users/nicco/Desktop/myPapers/ROAD_TO_THIRD_PAPER/The journals of gerontology/Rebuttal/final_files/figure_3.png", width = 5, height = 5, units = "in", res = 300)
cairo_ps(file = "/Users/nicco/Desktop/myPapers/ROAD_TO_THIRD_PAPER/The journals of gerontology/Rebuttal/final_files/figure_3.eps", onefile = F, fallback_resolution = 600, height = 7.5, width = 7.5)
layout(matrix(c(1,1, 1,1, 1,1, 2,3, 2,3), nrow = 5, ncol = 2, byrow = T))
revigoPlot(d)
lolliPlot_snp(traits, all.m)
lolliPlot_gene(traits_genes, overlapping.genes, geneList)
dev.off()
