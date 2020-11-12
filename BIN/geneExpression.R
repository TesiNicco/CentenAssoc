#take rna expression old vs. young
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO
gset <- getGEO(GEO="GSE11882", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
#thr: >80 and <50
gsms <- paste0("XX11111X0XXX1XXX1X1X0XXXXX1XXXXXXXXXXX1XXX0XXXXX1X",
               "XX0XXXXXXX1XXX0XXX1XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "1XXX0XXX0XXX0XXXXXXXXX0XXX0XXX0XX0XXX0XXX0XXX0XXX0",
               "XXX0XXXXXXX1XXXXXX1XXXX")
#thr: >80 and <30-65
gsms <- paste0("XX11111X0XXX1XXX1X1X0XXXXX1XXXXXXXXXXX1XXXXXXXXX1X",
               "XX0XXXXXXX1XXX0XXX1XXXXXXX0XXXXXXXXXXXXXXXXX0XXXXX",
               "1XXXXXXXXXXXXXXXXXXXXX0XXX0XXX0XXXXXX0XXXXXXX0XXX0",
               "XXX0XXXXXXX1XXXXXX1XXXX")

sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number="inf")

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
tT$newGene <- sapply(strsplit(tT$Gene.symbol, "///"), '[', 1)

#fromo here the new part on longevity genes -- 20191120
ind.info <- as.data.frame(gset$source_name_ch1)
xxx <- as.data.frame(str_split_fixed(ind.info$`gset$source_name_ch1`, ", ", 4))
ages <- str_split_fixed(xxx$V4, " ", 2)
xxx$V4 <- as.numeric(ages[, 1])
colnames(xxx) <- c("tissue", "region", "sex", "age")
xxx$stat <- 0
xxx$stat[which(xxx$age >= 80)] <- 1
table(xxx$stat)
mean(xxx$age[which(xxx$stat == 0)])
sd(xxx$age[which(xxx$stat == 0)])
mean(xxx$age[which(xxx$stat == 1)])
sd(xxx$age[which(xxx$stat == 1)])

#read longevity genes
genes.reported <- read.table("RESULTS/RESULTS_AnnotateMe_PrevRep/snp_annotation_geneList.txt", h=F)
genes.prs <- read.table("RESULTS/RESULTS_annotateMe_PRS/snp_annotation_geneList.txt", h=F)

# flag those reported and from the prs
genes.prs$source <- "PRS"
genes.reported$source <- "Reported"
all.genes <- data.frame(gene=c(as.character(genes.prs$V1), as.character(genes.reported$V1)), source=c(genes.prs$source, genes.reported$source))

# check which are in the expression data
info <- tT[, c("P.Value", "logFC", "newGene")]
sb <- merge(info, all.genes, by.x="newGene", by.y="gene")
sb <- sb[!duplicated(sb$newGene),]
sb <- sb[order(sb$P.Value),]
sb$adjP_mine <- p.adjust(sb$P.Value, method = "fdr", n=nrow(sb))
sb <- sb[order(sb$adjP_mine),]

# volcano plot
cairo_ps(file = "/Users/nicco/Desktop/myPapers/ROAD_TO_THIRD_PAPER/The journals of gerontology/Rebuttal/final_files/figure_4.eps", onefile = F, fallback_resolution = 600, height = 6.5, width = 6.5)
par(mfrow=c(1,1), mar=c(4, 4, 1, 5))

# separate prs from the reported based on pch -- 16:reported only -- 17:prs only -- 18:in common
sb$pch <- 16
#sb$pch[which(sb$source == "PRS")] <- 16
sb$pch[which(sb$newGene %in% genes.prs$V1 & sb$newGene %in% genes.reported$V1)] <- 18

# empty plot 
plot(sb$logFC, -log10(sb$P.Value), xlim=c(-2.5, 2.5), ylim=c(0, 8), pch=16, col="white", xlab="log2(Fold-change)",
     ylab="-log10(P-value)", cex.axis=0.75, cex.lab=1)

# grid
for (i in seq(-2.4, 2.4, 0.48)){abline(v=i, lwd=0.25, col="grey80")}
for (i in seq(0, 8, 0.8)){abline(h=i, lwd=0.25, col="grey80")}
abline(v=0, col="black", lwd=1)

# manage points size
sb$pSize <- 1
sb$pSize <- ceiling(-log10(sb$P.Value))
sb$pSize[which(sb$pSize >= 3)] <- 3
#colz <- colorRampPalette(c("navy", "blue", "light green", "dark green"))
colorz <- viridis(nrow(sb))
sb$col <- rev(colorz)

# put adj.pvalue line
tmp <- tail(sb[which(sb$adjP_mine <= 0.05), ], 1)
tmp.logP <- -log10(tmp$P.Value)
abline(h=tmp.logP, lwd=1, col="red", lty=2)

# then add points -- all for now
points(sb$logFC, -log10(sb$P.Value), pch=sb$pch, cex=sb$pSize, col=sb$col)
# then only those in common between prs and reported
points(sb$logFC[which(sb$pch==18)], -log10(sb$P.Value[which(sb$pch==18)]), pch=sb$pch[which(sb$pch==18)], cex=sb$pSize[which(sb$pch==18)], col="coral")

# text on top
text(x = 1.25, y = 8, labels = "Up in Old", adj = 0.5, cex=1.50, col="grey40", font=1)
text(x = -1.25, y = 8, labels = "Up in Young", adj = 0.5, cex=1.50, col="grey40", font=1)

# need to annotate
# define which to highlight -- top 10
highligh <- head(as.character(sb$newGene), 18)
overl <- as.character(head(sb[which(sb$pch == 18), "newGene"], 3))
highligh <- c(highligh, overl)
highligh <- highligh[!duplicated(highligh)]
pos <- c("ul", "ul", "ul", "r", "ul", "l", "r", "l", "ur", "d", "l", "ul", "dl", "r", "l", "dl", "r", "r")
x <- 0.4
for (i in 1:18){
  qq <- sb[which(sb$newGene == highligh[i]),]
  if (pos[i] == "u"){ text(x = qq$logFC, -log10(qq$P.Value)+x, labels = qq$newGene, font=3, cex=0.66, xpd=T); segments(x0 = qq$logFC, y0 = -log10(qq$P.Value)+x/4*3, x1 = qq$logFC, y1 = -log10(qq$P.Value), lwd=0.5) }
  if (pos[i] == "d"){ text(x = qq$logFC, -log10(qq$P.Value)-x/2, labels = qq$newGene, font=3, cex=0.66, xpd=T); segments(x0 = qq$logFC, y0 = -log10(qq$P.Value)-x/3, x1 = qq$logFC, y1 = -log10(qq$P.Value), lwd=0.5) }
  if (pos[i] == "ul"){ text(x = qq$logFC-x/2, -log10(qq$P.Value)+x, labels = qq$newGene, font=3, cex=0.66, xpd=T); segments(x0 = qq$logFC-x/2, y0 = -log10(qq$P.Value)+x/4*3, x1 = qq$logFC, y1 = -log10(qq$P.Value), lwd=0.5) }
  if (pos[i] == "dl"){ text(x = qq$logFC-x/2, -log10(qq$P.Value)-x, labels = qq$newGene, font=3, cex=0.66, xpd=T); segments(x0 = qq$logFC-x/2, y0 = -log10(qq$P.Value)-x/4*3, x1 = qq$logFC, y1 = -log10(qq$P.Value), lwd=0.5) }
  if (pos[i] == "ur"){ text(x = qq$logFC+x/2, -log10(qq$P.Value)+x, labels = qq$newGene, font=3, cex=0.66, xpd=T); segments(x0 = qq$logFC+x/2, y0 = -log10(qq$P.Value)+x/4*3, x1 = qq$logFC, y1 = -log10(qq$P.Value), lwd=0.5) }
  if (pos[i] == "uurr"){ text(x = qq$logFC+x, -log10(qq$P.Value)+x*1.75, labels = qq$newGene, font=3, cex=0.75, xpd=T); segments(x0 = qq$logFC+x, y0 = -log10(qq$P.Value)+x/2*3, x1 = qq$logFC, y1 = -log10(qq$P.Value), lwd=0.5) }
  if (pos[i] == "l"){ text(x = qq$logFC-x/4*3, -log10(qq$P.Value), labels = qq$newGene, font=3, cex=0.66, xpd=T); segments(x0 = qq$logFC-x/3, y0 = -log10(qq$P.Value), x1 = qq$logFC, y1 = -log10(qq$P.Value), lwd=0.5) }
  if (pos[i] == "r"){ text(x = qq$logFC+x, -log10(qq$P.Value), labels = qq$newGene, font=3, cex=0.66, xpd=T); segments(x0 = qq$logFC+x/3, y0 = -log10(qq$P.Value), x1 = qq$logFC, y1 = -log10(qq$P.Value), lwd=0.5) }
}
# legend
text(x = 3.1, y = 8 , labels = "-Log10(P)", cex=0.66, xpd=T, font=2)
gradient.rect(xleft = 3-0.125, ybottom = 6.75, xright = 3+0.125, ytop = 7.75, col=colorz, nslices = length(colorz), gradient = "horizontal", border = NA)
text(x = 3.15, y = 7.7, labels = round(max(-log10(sb$P.Value)), 1), adj = 0, xpd=T, font=1, cex=0.66)
text(x = 3.15, y = 6.85, labels = round(min(-log10(sb$P.Value)), 1), adj = 0, xpd=T, font=1, cex=0.66)
text(x = 3.15, y = 7.275, labels = round(median(-log10(sb$P.Value)), 1), adj = 0, xpd=T, font=1, cex=0.66)
points(x = 3, y = 5.8, pch=18, cex=2, col="coral", xpd=T)
text(x = 3, y = 6.2, labels = "Known\ngenes", cex=0.66, xpd=T, font=1)
dev.off()  




