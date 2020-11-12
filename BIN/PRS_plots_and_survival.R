###########################################
# THIS WILL PRODUCE PLOTS AND PERFORM SURVIVAL ANALYSIS ON THE BEST MODEL WITHOUT APOE
# LIBRARIES
library(stringr)
library(ggsci)
library(RColorBrewer)
library(data.table)
library(stringr)
library(survival)
library(survminer)

# PLOT SIGNIFICANCE OF PRS ACROSS THRESHOLDS FOR VARIANT INCLUSION
## read input
d <- read.table("RESULTS/PRS/all_models_stats_maf1pc.txt", h=T, dec=",")

bitmap("/Users/nicco/Desktop/myPapers/ROAD_TO_THIRD_PAPER/The journals of gerontology/Rebuttal/final_files/figure_1.tiff", width = 5, height = 3.5, units = "in", res = 300, type ="tifflzw")
postscript("/Users/nicco/Desktop/myPapers/ROAD_TO_THIRD_PAPER/The journals of gerontology/Rebuttal/final_files/figure_1.eps", width = 5, height = 4, horizontal = F, paper = "special", onefile = F)
colors = pal_npg("nrc", alpha = 0.7)(9)
colors <- brewer.pal(n = 9, name = "Set1")
colors_dots <- rep("grey80", 9)
colors_bars <- rep("deepskyblue3", 9)
par(mar=c(1, 4, 1, 5))
layout(matrix(c(1, 1, 1,
                2, 2, 2,
                2, 2, 2), nrow=3, byrow=TRUE))

plot(0, 0, pch=16, col="white", xlim=c(1, 9), ylim=c(0, 5), xaxt="none", 
     yaxt='none', xlab="", ylab=expression("# of SNPs"), cex.lab=1, bty='n', xpd=T)
for (x in seq(1, 10, 10/10)){abline(v=x, lwd=0.5, col="grey80")}
for (x in seq(0, 5, 5/4)){abline(h=x, lwd=0.5, col="grey80")}
for (x in 1:nrow(d)){
  rect(xleft = x-0.20, ybottom = 0, xright = x+0.20, ytop = log10(d$number_variants[x]), 
       col = colors_bars[x], lwd=0.60)
  text(x = x, y = log10(d$number_variants[x])+0.5, cex=0.66, labels = d$number_variants[x], xpd=T)
  }
par(mar=c(4, 4, 0, 5))
plot(0, 0, pch=16, col="white", xlim=c(1, 9), ylim=c(0, 12), xaxt="none", 
     xlab="PRS", ylab="-log10(P-value)", cex.lab=1, bty='n', cex.axis=0.75)
for (x in seq(1, 10, 10/10)){abline(v=x, lwd=0.25, col="grey80")}
for (x in seq(0, 12, 12/10)){abline(h=x, lwd=0.25, col="grey80")}
labs <- c("Known", "PRS-8", "PRS-7", "PRS-6", "PRS-5", "PRS-4", "PRS-3", "PRS-2", "PRS-1")
for (i in 1:9){text(x = i, y = -1, labels = labs[i], adj = 1, xpd=T, cex=0.75, srt=45)}
points(seq(1, 9), -log10(d$P_apoe), col=colors_dots, pch=16, cex=1.25, type = "b") 
points(seq(1, 9), -log10(d$P), col=colors_dots, pch=17, cex=1.25, type="b")
# put the interesting dots in red
points(x=5, y=-log10(d$P[5]), col="red", pch=17, cex=1.25, type = "b") 
points(x=4, y=-log10(d$P_apoe[4]), col="red", pch=16, cex=1.25, type="b")
legend(x = 9.25, y = 12, legend = c("APOE Inc.", "APOE Exc."), 
       pch=c(16, 17), col="black", cex=0.80, pt.cex = 1, xpd=T, bty="n")
dev.off()

# PLOT ROC CURVE FOR THE DIFFERENT PRS
## read roc object
load("RESULTS/PRS/ROC_curve_objects.RData")
roc <- all.roc[[1]]
roc.noAPOE <- all.roc[[2]]

function.plotROC_pPRS_PRS <- function(l, l.noAPOE){
  #graphical parameters
  par(mfrow=c(2, 1), mar=c(4, 5, 4, 8))
  
  #color palette
  colors <- brewer.pal(n = 9, name = "Set1")
  
  #AD plot with APOE
  #empty plot as background
  plot(x = 0, y = 0, type = "p", pch=16, col="white", xlim = c(1, 0), ylim=c(0, 1), 
       xlab="Specificity", ylab="Sensitivity", cex.lab=1.50, main="PRS ~ APOE included",
       cex.axis=1.25, cex.main=2)
  
  #some grid
  for (i in seq(0, 1, 0.1)){abline(h=i, lty=1, lwd=0.5, col="grey80"); abline(v=i, lty=1, lwd=0.5, col="grey80")}
  
  #diagunal line for reference
  segments(x0 = 1, y0 = 0, x1 = 0, y1 = 1, col = "red", lwd=1.25)
  
  #then loop for the pPRS aucs
  for (i in 1:length(l)){
    lines(l[[i]], col=colors[i], lty=1, lwd=2)
  }
  legend(x = -0.1, y = 1.05, legend = c("Reported", "PRS8", "PRS7", "PRS6", "PRS5", "PRS4", "PRS3", "PRS2", "PRS1"), lty=1, 
         col=c(colors[1], colors[2], colors[3], colors[4], colors[5], colors[6], colors[7], colors[8], colors[9]), cex=1.10, xpd=T, bty="n", lwd=3)
  
  #AD plot without APOE
  par(mar=c(4, 5, 4, 8))
  #empty plot as background
  plot(x = 0, y = 0, type = "p", pch=16, col="white", xlim = c(1, 0), ylim=c(0, 1), 
       xlab="Specificity", ylab="Sensitivity", cex.lab=1.50, main="PRS ~ APOE excluded",
       cex.axis=1.25, cex.main=2)
  
  #some grid
  for (i in seq(0, 1, 0.1)){abline(h=i, lty=1, lwd=0.5, col="grey80"); abline(v=i, lty=1, lwd=0.5, col="grey80")}
  
  #diagunal line for reference
  segments(x0 = 1, y0 = 0, x1 = 0, y1 = 1, col = "red", lwd=1.25)
  
  #then loop for the pPRS aucs
  for (i in 1:length(l.noAPOE)){
    lines(l.noAPOE[[i]], col=colors[i], lty=1, lwd=2)
  }
  
  #now legend
  legend(x = -0.1, y = 1.05, legend = c("PRS8", "PRS7", "PRS6", "PRS5", "PRS4", "PRS3", "PRS2", "PRS1"), lty=1, 
         col=c(colors[1], colors[2], colors[3], colors[4], colors[5], colors[6], colors[7], colors[8]), cex=1.10, xpd=T, bty="n", lwd=3)
}
library(ggplot2)
pdf("RESULTS/PRS/ROC_PRS_maf1pc_v1.pdf", height = 10, width = 7)
function.plotROC_pPRS_PRS(roc, roc.noAPOE)
dev.off()

# SURVIVAL ANALYSIS
## function to plot kaplan meier with left and right truncation
function.SurvivalCurve_left_right_truncated <- function(fit1, cox.prs){
  #get data for the plot
  #survival estimates
  survival.estimates <- fit1$surv
  #times
  times <- fit1$time + min(dat$age_baseline)
  #confidence intervals
  ci.lower <- fit1$lower
  ci.upper <- fit1$upper
  #numbers at risk
  n.risk <- fit1$n.risk
  
  #put data together
  #need to find the splitting point somewhere in the middle of the data
  flipp.point <- which(survival.estimates == 1)
  flipp.point <- (flipp.point[which(flipp.point > 200)][1])-1
  gr1 <- data.frame(surv = survival.estimates[1:flipp.point], time = times[1:flipp.point],
                    lower = ci.lower[1:flipp.point], upper = ci.upper[1:flipp.point], 
                    n = n.risk[1:flipp.point])
  start <- flipp.point+1
  gr2 <- data.frame(surv = survival.estimates[start:length(survival.estimates)], 
                    time = times[start:length(survival.estimates)],
                    lower = ci.lower[start:length(survival.estimates)],
                    upper = ci.upper[start:length(survival.estimates)],
                    n = n.risk[start:length(survival.estimates)])
  
  #HERE IS THE PLOT
  #graphical parameters
  layout(matrix(c(1, 1, 1, 2), nrow=4))
  par(mar=c(5, 8, 2, 4))
  
  #select minimum and maximum
  gl.min <- min(gr1$time, gr2$time)
  gl.max <- max(gr1$time, gr2$time)
  gl.min <- 55
  gl.max <- 105
  
  #empty plot
  plot(0, 0, type="s", col="white", cex.lab=2, cex.axis=1.50, xlab="Age",
       ylab="Survival probabilities", xlim=c(gl.min, gl.max), ylim=c(0, 1))
  
  #grid
  for (i in seq(gl.min, gl.max, (gl.max-gl.min)/10)){abline(v=i, lwd=0.5, col="grey80")}
  for (i in seq(0, 1, 0.10)){abline(h=i, lwd=0.5, col="grey80")}
  
  #then the curves: use polygon to color confidence interval: first create the vectors for the polygon
  #make sure to exclude NA
  x <- c(gr1$time, rev(gr1$time))
  y <- c(gr1$lower, rev(gr1$upper))
  y[which(is.na(y))] <- 0
  x2 <- c(gr2$time, rev(gr2$time))
  y2 <- c(gr2$lower, rev(gr2$upper))
  y2[which(is.na(y2))] <- 0
  
  #then draw polygon
  polygon(x, y, col = alpha("yellow", 0.5), border = NA)
  polygon(x2, y2, col = alpha("deepskyblue2", 0.5), border = NA)
  
  #then put lines
  lines(gr1$time, gr1$surv, type="s", lwd=3, col="gold3")
  lines(gr2$time, gr2$surv, type="s", lwd=3, col="navy")
  
  #need to draw dotted line at p=0.5: for this need to find intersection points
  gr1$diff <- abs(0.5-gr1$surv)
  gr1 <- gr1[order(gr1$diff),]
  gr2$diff <- abs(0.5-gr2$surv)
  gr2 <- gr2[order(gr2$diff),]
  segments(x0 = 0, y0 = 0.5, x1 = gr2$time[1], y1 = 0.5, lwd=1.5, lty=2, col="black")
  segments(x0 = gr1$time[1], y0 = 0.5, x1 = gr1$time[1], y1 = 0, lwd=1.5, lty=2, col="black")
  segments(x0 = gr2$time[1], y0 = 0.5, x1 = gr2$time[1], y1 = 0, lwd=1.5, lty=2, col="black")
  
  #pvalue on bottom-right
  text(x = 55, y = 0.1, labels = paste(expression(p), " = ", round(summary(cox.prs)$coefficients[1, 5], 3), sep=""), cex=2.5, adj=0, font=3)
  
  #legend on top-right
  legend("topright", legend = c("Low-PRS", "High-PRS"), lwd=2, col=c("gold3", "navy"), bty="n", cex=2)
  
  #second plot: numbers at risk
  par(mar=c(2, 8, 2, 4))
  
  #empty plot
  plot(0, 0, type="s", col="white", cex.lab=1.75, cex.axis=1.50, xlab="",
       ylab="", xlim=c(gl.min, gl.max), ylim=c(0, 0.5), xaxt="n", yaxt='n')
  
  #grid
  for (i in seq(gl.min, gl.max, (gl.max-gl.min)/10)){abline(v=i, lwd=0.5, col="grey80")}
  for (i in seq(0, 1, 0.25)){abline(h=i, lwd=0.5, col="grey80")}
  
  #choose where to put numbers
  where <- c(60, 70, 80, 90, 100)
  gr1 <- gr1[order(gr1$time),]
  gr2 <- gr2[order(gr2$time),]
  for (i in where){
    #get data
    n1.all <- gr1$n[which(gr1$time >= i)]
    n1 <- n1.all[1]
    n2.all <- gr2$n[which(gr2$time >= i)]
    n2 <- n2.all[1]
    
    #put data in
    text(x = i, y = 0.125, labels = n2, font=2, cex=2)
    text(x = i, y = 0.375, labels = n1, font=2, cex=2)
  }
  
  #finally put labels on the left
  text(x = 52.5, y = 0.125, labels = "High-PRS", col="navy", font=2, adj=1, xpd=T, cex=1.75)
  text(x = 52.5, y = 0.375, labels = "Low-PRS", col="gold3", font=2, adj=1, xpd=T, cex=1.75)
  text(x = 55, y = 0.55, labels = "Numbers at risk", col="black", cex=1.75, adj=0, xpd=T, font=2)  
  #dev.off()
}

## same function as the previous but with stratification by apoe
function.SurvivalCurve_left_right_truncated_APOE <- function(fit1, cox.binary.apoeYesNo){
  #get data for the plot
  #survival estimates
  survival.estimates <- fit1$surv
  #times
  times <- fit1$time + min(dat$age_baseline)
  #confidence intervals
  ci.lower <- fit1$lower
  ci.upper <- fit1$upper
  #numbers at risk
  n.risk <- fit1$n.risk
  
  #put data together
  #need to find the splitting point somewhere in the middle of the data
  flipp.point <- which(survival.estimates == 1)
  #identify breaking points
  stop = c()
  i <- 1
  xx <- flipp.point[1]
  while (i <= (length(flipp.point)-1)){
    if (flipp.point[i+1] == xx+1){
      i <- i +1
      xx <- xx +1
    } else {
      stop <- c(stop, flipp.point[i+1])
      i <- i+1
      xx <- flipp.point[i]
    }
  }
  stop <- c(stop, length(survival.estimates))
  
  #output has to be a list
  l = list()
  #manually input flipping points
  for (i in 1:length(stop)){
    v <- stop[i]
    if (i == 1){ 
      start <- 1
      end <- v-1
    } else if (i == length(stop)){
      start <- stop[i-1]
      end <- v
    } else {
      start <- stop[i-1]
      end <- v-1 
    }
    
    #create dataframe
    tmp <- data.frame(surv = survival.estimates[start:end], time = times[start:end],
                              lower = ci.lower[start:end], upper = ci.upper[start:end], 
                              n = n.risk[start:end])
    #add to list
    l[[i]] <- tmp
    
  }
  
  #HERE IS THE PLOT
  #graphical parameters
  layout(matrix(c(1,1, 2), nrow=3))
  par(mar=c(4, 5, 2, 1))
  
  #select minimum and maximum
  gl.min <- min(fit1$time)
  gl.max <- max(fit1$surv)
  gl.min <- 55
  gl.max <- 105
  
  #empty plot
  plot(0, 0, type="s", col="white", cex.lab=1, cex.axis=0.75, xlab="Age",
       ylab="Survival probabilities", xlim=c(gl.min, gl.max), ylim=c(0, 1), bty="n")
  
  #grid
  for (i in seq(gl.min, gl.max, (gl.max-gl.min)/10)){abline(v=i, lwd=0.25, col="grey80")}
  for (i in seq(0, 1, 0.10)){abline(h=i, lwd=0.25, col="grey80")}
  
  #then draw polygon
  library(ggsci)
  colz <- pal_lancet(palette = "lanonc")(9)
  colz <- brewer.pal(9, "Blues")
  colz <- colz[c(3, 5, 7, 9)]
  for (i in 1:length(l)){ 
    df <- l[[i]]
    points(df$time, df$surv, col = colz[i], type="s", lwd=1) 
    }

  low.nonCarr <- l[[1]]
  low.Carr <- l[[2]]
  high.nonCarr <- l[[3]]
  high.Carr <- l[[4]]
  
  #need to draw dotted line at p=0.5: for this need to find intersection points
  low.nonCarr$diff <- abs(0.5-low.nonCarr$surv)
  low.nonCarr <- low.nonCarr[order(low.nonCarr$diff),]
  low.Carr$diff <- abs(0.5-low.Carr$surv)
  low.Carr <- low.Carr[order(low.Carr$diff),]
  high.nonCarr$diff <- abs(0.5-high.nonCarr$surv)
  high.nonCarr <- high.nonCarr[order(high.nonCarr$diff),]
  high.Carr$diff <- abs(0.5-high.Carr$surv)
  high.Carr <- high.Carr[order(high.Carr$diff),]
  segments(x0 = 0, y0 = 0.5, x1 = high.nonCarr$time[1], y1 = 0.5, lwd=0.75, lty=2, col="black")
  segments(x0 = low.nonCarr$time[1], y0 = 0.5, x1 = low.nonCarr$time[1], y1 = 0, lwd=0.75, lty=2, col="black")
  segments(x0 = low.Carr$time[1], y0 = 0.5, x1 = low.Carr$time[1], y1 = 0, lwd=0.75, lty=2, col="black")
  segments(x0 = high.Carr$time[1], y0 = 0.5, x1 = high.Carr$time[1], y1 = 0, lwd=0.75, lty=2, col="black")
  segments(x0 = high.nonCarr$time[1], y0 = 0.5, x1 = high.nonCarr$time[1], y1 = 0, lwd=0.75, lty=2, col="black")
  
  # put p50 into a dataframe
  p50_ages <- data.frame(Group=c("low_Carr", "low_nonCarr", "high_Carr", "high_nonCarr"), age=c(low.Carr$time[1], low.nonCarr$time[1], high.Carr$time[1], high.nonCarr$time[1]))
  list_p50 <- list(low.nonCarr$time[1], low.Carr$time[1], high.Carr$time[1], high.nonCarr$time[1])
  #pvalue on bottom-right
  text(x = 55, y = 0.1, labels = paste(expression(p), " = ", round(summary(cox.binary.apoeYesNo)$coefficients[1, 5], 3), sep=""), cex=1, adj=0, font=3)
  
  #legend on top-right
  legend("topright", legend = c("Low-PRS", "Low-PRS APOE4+", "High-PRS", "High-PRS APOE4+"), lwd=1, col=colz[1:4], bty="n", cex=0.75)
  
  #second plot: numbers at risk
  par(mar=c(2, 5, 1, 1))
  
  #empty plot
  plot(0, 0, type="s", col="white", cex.lab=1.75, cex.axis=1.50, xlab="",
       ylab="", xlim=c(gl.min, gl.max), ylim=c(0, 0.5), xaxt="n", yaxt='n', bty="n")
  
  #grid
  for (i in seq(gl.min, gl.max, (gl.max-gl.min)/10)){abline(v=i, lwd=0.25, col="grey80")}
  for (i in seq(0, 0.5, (0.5/4))){abline(h=i, lwd=0.25, col="grey80")}
  
  #choose where to put numbers
  where <- c(60, 70, 80, 90, 100)
  low.nonCarr <- low.nonCarr[order(low.nonCarr$time),]
  low.Carr <- low.Carr[order(low.Carr$time),]
  high.nonCarr <- high.nonCarr[order(high.nonCarr$time),]
  high.Carr <- high.Carr[order(high.Carr$time),]
  for (i in where){
    #get data
    n1.all <- low.nonCarr$n[which(low.nonCarr$time >= i)]
    n1 <- n1.all[1]
    n2.all <- low.Carr$n[which(low.Carr$time >= i)]
    n2 <- n2.all[1]
    n3.all <- high.nonCarr$n[which(high.nonCarr$time >= i)]
    n3 <- n3.all[1]
    n4.all <- high.Carr$n[which(high.Carr$time >= i)]
    n4 <- n4.all[1]
    
    #put data in
    text(x = i, y = 0.10, labels = n2, font=1, cex=0.80, col=colz[2])
    text(x = i, y = 0.20, labels = n1, font=1, cex=0.80, col=colz[1])
    text(x = i, y = 0.30, labels = n4, font=1, cex=0.80, col=colz[4])
    text(x = i, y = 0.40, labels = n3, font=1, cex=0.80, col=colz[3])
  }
  
  #finally put labels on the left
  text(x = 55, y = 0.51, labels = "Numbers at risk", col="black", cex=0.80, adj=0, xpd=T, font=1)  
  text(x = 52.5, y = 0.1, labels = "Low-PRS\nAPOE4+", col=colz[2], font=1, adj=1, xpd=T, cex=0.70)
  text(x = 52.5, y = 0.2, labels = "Low-PRS", col=colz[1], font=1, adj=1, xpd=T, cex=0.70)
  text(x = 52.5, y = 0.3, labels = "High-PRS\nAPOE4+", col=colz[4], font=1, adj=1, xpd=T, cex=0.70)
  text(x = 52.5, y = 0.4, labels = "High-PRS", col=colz[3], font=1, adj=1, xpd=T, cex=0.70)
  #dev.off()
}

# MAIN
## find best model
d <- d[order(d$P), ]
mostSign <- rownames(d)[1]
## get PRS of that model
fname <- paste("RESULTS/PRS/prs_", mostSign, "_noAPOE.txt", sep="")
prs <- fread(fname, h=T)

## read phenotypes
pheno <- read.table("INPUTS_OTHER/lasa_phenotypes_QCpassed.table", h=T, stringsAsFactors=F, sep="\t")
## reduce phenotype file
pheno$surv_time <- pheno$age_GWA - pheno$age_baseline
pheno <- pheno[, c("respnr", "death", "surv_time", "age_baseline", "sex")]
pheno <- pheno[which(pheno$death != "-1"),]
pheno$respnr <- as.character(pheno$respnr)
## merge phenotypes and prs
dat <- merge(prs, pheno, by.x="IID", by.y="respnr")
## create survival object: left-censored, right-truncated
surv_object <- Surv(time = dat$age_baseline-min(dat$age_baseline), time2 = dat$age_baseline+dat$surv_time-min(dat$age_baseline),
            event = dat$death, type="counting")
## split the population based on the PRSs
train <- dat[which(dat$age_baseline <= 65),]
## quantiles
q.n <- 2
quant.vec <- seq(0+1/q.n, 1-1/q.n, 1/q.n)
dat$Quant <- 1
counter <- 2
for (x in quant.vec){
    q <- quantile(train$PRS, probs = x)
    ## then assign to bigger dataset (m1)
    dat$Quant[which(dat$PRS >= q)] <- counter
    counter <- counter + 1
}
table(dat$Quant)

## exploration here -- try to see whether the HR is lower or larger for females
apoe.gt <- read.table("INPUTS_OTHER/apoe_genotype_everyone.txt", h=T)
# add apoe to genotype
dat <- merge(dat, apoe.gt, by="IID")
dat$APOE4 <- "no"
dat$APOE4[which(dat$apoe_geno %in% c("e4/e3", "e4/e4", "e4/e2"))] <- "yes"

## cox models
cox.prs <- coxph(surv_object ~ dat$PRS + dat$PC1 + dat$PC2 + dat$PC3 + dat$PC4 + dat$PC5 + dat$sex)
summary(cox.prs)

## survival object for kaplan meier  
fit1 <- survfit(surv_object ~ Quant, data = dat)  

## plots
function.SurvPlot <- function(dat, t1, t2, fit1, cox.prs, surv_object){
  #apoe.gt <- read.table("GENO_DATA/apoe_genotype_everyone.txt", h=T)
  
  # add apoe to genotype
  dat <- merge(dat, apoe.gt, by="IID")
  dat$APOE4 <- "no"
  dat$APOE4[which(dat$apoe_geno %in% c("e4/e3", "e4/e4", "e4/e2"))] <- "yes"
  
  # make the plot without APOE
  pdf(t1, height=10, width=10)
  function.SurvivalCurve_left_right_truncated(fit1, cox.prs)
  dev.off()
  
  # also the plot with apoe stratification
  cox.binary.apoeYesNo <- coxph(surv_object ~ dat$PRS + dat$PC1 + dat$PC2 + dat$PC3 + dat$PC4 + dat$PC5 + dat$sex + dat$APOE4)
  fit1 <- survfit(surv_object ~ Quant + APOE4, data = dat)
  
  # interaction term of prs and gender
  cox.binary.apoeYesNo.interaction <- coxph(surv_object ~ dat$PRS * dat$sex + dat$PC1 + dat$PC2 + dat$PC3 + dat$PC4 + dat$PC5 + dat$APOE4)
  summary(cox.binary.apoeYesNo.interaction)
  cox.binary.apoeYesNo.interaction <- coxph(surv_object ~ dat$PRS * dat$APOE4 + dat$PC1 + dat$PC2 + dat$PC3 + dat$PC4 + dat$PC5 + dat$sex)
  summary(cox.binary.apoeYesNo.interaction)
  
  pdf(t2, height=10, width=10)
  #pdf("RESULTS/PRS/survPlot_strat.pdf", height=4, width=3.4)
  #png("RESULTS/PRS/survPlot_strat.png", height=4, width=3.4, units = "in", res=1200)
  bitmap("/Users/nicco/Desktop/myPapers/ROAD_TO_THIRD_PAPER/The journals of gerontology/Rebuttal/final_files/figure_2.tiff", width = 3.5, height = 5, units = "in", res = 300, type ="tifflzw")
  postscript("/Users/nicco/Desktop/myPapers/ROAD_TO_THIRD_PAPER/The journals of gerontology/Rebuttal/final_files/figure_2.eps", width = 4, height = 5, horizontal = F, paper = "special", onefile = F)
  function.SurvivalCurve_left_right_truncated_APOE(fit1, cox.binary.apoeYesNo)
  dev.off()
}
function.SurvPlot(dat, "RESULTS/PRS/survPlot_noAPOE.pdf", "RESULTS/PRS/survPlot_APOE.pdf", fit1, cox.prs, surv_object)

# final cox
cox.binary.apoeYesNo <- coxph(surv_object ~ dat$PRS + dat$PC1 + dat$PC2 + dat$PC3 + dat$PC4 + dat$PC5 + dat$sex + dat$APOE4)
fit1 <- survfit(surv_object ~ Quant + APOE4, data = dat)
summary(cox.binary.apoeYesNo)

males <- dat[which(dat$sex == "male"), ]
females <- dat[which(dat$sex == "female"), ]
## make separate survival objects
surv_object_males <- Surv(time = males$age_baseline-min(males$age_baseline), time2 = males$age_baseline+males$surv_time-min(males$age_baseline),
            event = males$death, type="counting")
surv_object_females <- Surv(time = females$age_baseline-min(females$age_baseline), time2 = females$age_baseline+females$surv_time-min(females$age_baseline),
            event = females$death, type="counting")
## make separate cox models for males and females
cox.males <- coxph(surv_object_males ~ males$PRS + males$PC1 + males$PC2 + males$PC3 + males$PC4 + males$PC5 + males$APOE4)
summary(cox.males)
cox.females <- coxph(surv_object_females ~ females$PRS + females$PC1 + females$PC2 + females$PC3 + females$PC4 + females$PC5 + females$APOE4)
summary(cox.females)
## finally, for the plot, make a function to do a forestplot
forrestPlot <- function(cox.males, cox.females, cox.prs){
  # define x and y limits
  xlm <- c(-1, 3)
  ylm <- c(0, 4)
  # graphical settings
  par(mfrow=c(1, 1), mar=c(4, 5, 4, 5))
  # basic background plot
  plot(0, 0, pch=16, col="white", bty='none', xaxt='n', yaxt='n', xlab="", ylab="", xlim=xlm, ylim=ylm)
  # grid
  for (x in seq(0, 2, 0.5)){ abline(v=x, lwd=0.4, col="grey80") }
  for (x in seq(0, 4, 0.75)){ segments(x0 = 0, y0= x, x1 = 2, y1 = x, lwd=0.4, col="grey80") }
  # vertical line at 0
  abline(v=1, lwd=1.5, lty=1, col="black")
  text(x=1, y=4.4, labels="Hazard Ratio of PRS", cex=1.50, font=2, xpd=T)
  # put x-axis
  axis(side=1, at=seq(0, 2, 0.5), labels=seq(0, 2, 0.5), cex.axis=1.5)
  # put points
  segments(x0=exp(summary(cox.males)$coefficients[1, 1]-(1.96*summary(cox.males)$coefficients[1, 3])), y0=1, x1=exp(summary(cox.males)$coefficients[1, 1]+(1.96*summary(cox.males)$coefficients[1, 3])), y1=1, lwd=2.5, col="navy")
  points(x=exp(cox.males$coefficients[1]), y=1, pch=17, col="navy", cex=1.75)
  segments(x0=exp(summary(cox.females)$coefficients[1, 1]-(1.96*summary(cox.females)$coefficients[1, 3])), y0=2, x1=exp(summary(cox.females)$coefficients[1, 1]+(1.96*summary(cox.females)$coefficients[1, 3])), y1=2, lwd=2.5, col="pink")
  points(x=exp(cox.females$coefficients[1]), y=2, pch=18, col="pink", cex=2.5)  
  segments(x0=exp(summary(cox.prs)$coefficients[1, 1]-(1.96*summary(cox.prs)$coefficients[1, 3])), y0=3, x1=exp(summary(cox.prs)$coefficients[1, 1]+(1.96*summary(cox.prs)$coefficients[1, 3])), y1=3, lwd=2.5, col="coral")
  points(x=exp(cox.prs$coefficients[1]), y=3, pch=16, col="coral", cex=1.9)
  # some text
  text(x=-1.5, y=1, adj=0, labels=paste("Males (N=", cox.males$n, ")", sep=""), col="navy", cex=1.50, font=2, xpd=T)
  text(x=2.5, y=1, adj=0, labels=paste("p=", round(summary(cox.males)$coefficients[1, 5], 2), sep=""), col="navy", cex=1.20, font=3, xpd=T) 
  text(x=-1.5, y=2, adj=0, labels=paste("Females (N=", cox.females$n, ")", sep=""), col="pink", cex=1.50, font=2, xpd=T)
  text(x=2.5, y=2, adj=0, labels=paste("p=", round(summary(cox.females)$coefficients[1, 5], 2), sep=""), col="pink", cex=1.20, font=3, xpd=T) 
  text(x=-1.5, y=3, adj=0, labels=paste("Combined (N=", cox.males$n+cox.females$n, ")", sep=""), col="coral", cex=1.50, font=2, xpd=T)
  text(x=2.5, y=3, adj=0, labels=paste("p=", round(summary(cox.prs)$coefficients[1, 5], 2), sep=""), col="coral", cex=1.20, font=3, xpd=T) 
  
}
pdf("RESULTS/PRS_ANNOTATION/forrestPlot_male_females_survival.pdf", height=4, width=8)
forrestPlot(cox.males, cox.females, cox.binary.apoeYesNo)
dev.off()

############################################
############################################
# compare hr of PRS in males and females -- this is done in function: function.SurvPlot
