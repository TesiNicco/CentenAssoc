###################################################
# SCRIPT TO PERFORM PRS AT DIFFERENT PVALUE THRESHOLD DEFINED IN THE CODE
# LIBRARIES
library(stringr)
library(ggplot2)
library(ggsci)
library(data.table)
library(parallel)
library(lme4)
library(pROC)
args = commandArgs(trailingOnly=TRUE)

# FUNCTIONS
## function to calculate PRS and also do the association of single variants
function.PRS_mod <- function(dos.ph, sb, assoc){
        prs <- as.data.frame(matrix(data=0, nrow=nrow(dos.ph), ncol=3))
        colnames(prs) <- c("IID", "PRS", "PHENO")
        prs$IID <- dos.ph[, "IID"]
        prs$PHENO <- dos.ph[, "PHENO"]
        n.var <- 0
        Max <- ncol(dos.ph)-2
        snps.assoc <- as.data.frame(matrix(data=NA, nrow=(Max-2), ncol=6))
        colnames(snps.assoc) <- c("snp", "a1", "beta", "se", "pvalue", "direction")
        covar <- read.table("GENO_DATA/20200123_covariates_kept_Age_gwas.txt", h=T, stringsAsFactors=F)
        
        samp.names <- dos.ph$IID
        dos.ph$IID <- NULL
        dos.ph <- as.matrix(dos.ph)

        for (snp in seq(1, Max)){

                #select snp
                snp.data <- dos.ph[, snp]

                #snp name
                snp.name <- colnames(dos.ph)[snp]
                snp.name2 <- str_split_fixed(snp.name, "_", 2)
                loc = snp.name2[, 1]
                a1 <- snp.name2[, 2]

                #check beta
                if (length(grep(":", loc)) > 0){
                        check <- sb[which(sb$LOCUS == loc),]
                } else {
                        check <- sb[which(sb$SNP == loc),]
                }

                if (nrow(check) != 0){
                    n.var <- n.var + 1
                    if (a1 != toupper(check$A1)){ snp.data <- 2 - snp.data; a1 <- check$A1}
                    if (check$BETA < 0){ check$BETA <- check$BETA * (-1); snp.data <- 2 - snp.data; a1 <- toupper(check$A2) }
                    if (check$BETA < 0){ print("!!! Problem !!!") }

                    #get tmp score
                    tmp <- check$BETA * snp.data
                    prs$PRS <- prs$PRS + tmp
                }

                #if requested, also do association
                if (assoc =="yes"){
                    df <- data.frame(IID = samp.names, dos = snp.data, pheno = dos.ph[, "PHENO"]-1)
                    df <- merge(df, covar, by="IID")
                    m1 <- glm(pheno ~ dos + PC1 + PC2 + PC3 + PC4 + PC5, data=df, family="binomial")
                    snps.assoc[(snp), ] <- c(loc, a1, summary(m1)$coefficients[2, c(1, 2, 4)], NA)
                
                    #finally, check positions
                    if ((as.numeric(snps.assoc$beta[snp]) * check$BETA) >= 0) { snps.assoc[snp, "direction"] <- "same" } else { snps.assoc[snp, "direction"] <- "opposite"}
                }
        }
        prs$IID <- unlist(prs$IID)
        prs$PHENO <- unlist(prs$PHENO)
        print(paste("    PRS done with ", n.var, sep=""))
        l <- list(prs, snps.assoc)
        return(l)
}

## function to extract dosages and prepare for the prs
function.extractAndPrepare_mod <- function(ukb.sbs, maf_thr){
        my.match.rsid <- ukb.sbs
        #write temporary file
        write.table(my.match.rsid$SNP, "RESULTS/PRS/tmp.variants", quote=F, row.names=F, col.names=F)

        #run plink2 to extract dosages for the PRS
        chrom_list <- as.character(unique(my.match.rsid$CHR))
        chroms <- paste0(chrom_list, collapse=" ")
        cmd_plink <- paste("for chr in ", chroms, "; do plink2 --pfile GENO_DATA/PLINK2/chr${chr}_plink2_upd --extract RESULTS/PRS/tmp.variants --export A --out RESULTS/PRS/chr${chr}_dosages; done", sep="")
        system(cmd_plink, ignore.stdout = TRUE)
        #clean log files
        cmd_clean <- "rm RESULTS/PRS/*log"
        system(cmd_clean)
        #merge files together
        cmd_bind <- "paste --delimiters='\t' RESULTS/PRS/*raw > RESULTS/PRS/chrAll_dosages.txt"
        system(cmd_bind)
        cmd_clean <- "rm RESULTS/PRS/*raw"
        system(cmd_clean)
        cmd_clean <- "rm RESULTS/PRS/tmp*"
        system(cmd_clean)

        #read file back in R
        dos <- fread("RESULTS/PRS/chrAll_dosages.txt", h=T, sep="\t", check.names=F, stringsAsFactors=F)

        cmd_clean <- "rm RESULTS/PRS/chrAll_dosages.txt"
        system(cmd_clean)

        #take informative columns
        x <- colnames(dos)
        col.index <- grep(":", x)
        col.index2 <- grep("rs", x)
        all.index <- c(col.index, col.index2)
        dos.cl <- dos[, ..all.index]
        dos.cl$IID <- dos$IID
        dos.cl$PHENO <- dos$PHENOTYPE

        return(dos.cl)
}

# MAIN
## set maf threshold
maf_thr <- 0.01
do_assoc <- args[1]
#help message in case
if (do_assoc == "--help"){ print("Choices are: yes/no") }

# read summary statistics
print("## Reading Summary Statistics Data")
ukb <- fread("INPUTS_OTHER/GWAS_data_Timmers_cleaned.txt", h=T)      #with new pipeline

# read clumped variants
clumped <- fread("RESULTS/CLUMPING/clumping_chrAll_variantsKept.txt", h=T, stringsAsFactors=F)

# read reported variants
reported <- fread("RESULTS/SINGLE_VARIANT/reported_variants_info.txt", h=F, stringsAsFactors=F)
colnames(reported) <- c("CHR", "SNP", "BO", "BP", "REF", "ALT")

#what thresholds to make the PRS? --> 0==reported -- 5e-8 -- 5e-7 -- 5e-6 -- 5e-5 -- 5e-4 -- 5e-3 -- 5e-2 -- 5e-1
thrs <- c(0, 5e-8, 5e-7, 5e-6, 5e-5, 5e-4, 5e-3, 5e-2, 5e-1)

#general output with model statistics
model.stats <- as.data.frame(matrix(data=NA, nrow=length(thrs), ncol=7))
colnames(model.stats) <- c("number_variants", "OR_apoe", "SE_apoe", "P_apoe", "OR", "SE", "P")
rownames(model.stats) <- thrs
c <- 1
roc.list <- list()
roc.list.noAPOE <- list()
assoc.list <- list()

# main loop
for (thr in thrs){
        print(paste("### PRS with variant inclusion of ", thr, sep=""))

        #check if is the reported turn
        if (thr == 0){
            ukb.sbs <- ukb[which(ukb$SNP %in% reported$SNP), ]
            ukb.sbs$LOCUS <- paste(ukb.sbs$CHR, ukb.sbs$BP, sep=":")
            print(paste("    PRS should include ", nrow(ukb.sbs), " variants", sep=""))
        } else {
            #make subset on clumped variants
            clumped.sbs <- clumped[which(clumped$P <= thr),]
            tmp.apoe2 <- data.frame(CHR=19, SNP="rs7412", BP=1234, P=5e-100, TOTAL=1000)
            clumped.sbs <- rbind(clumped.sbs, tmp.apoe2)

            #take the same variants from the sumstats -- remember to add apoe2 as it is excluded normally
            ukb.sbs <- ukb[which(ukb$SNP %in% clumped.sbs$SNP),]
            ukb.sbs$LOCUS <- paste(ukb.sbs$CHR, ukb.sbs$BP, sep=":")
            print(paste("    PRS should include ", nrow(ukb.sbs), " variants", sep=""))
        }

        out_name <- paste("RESULTS/PRS/variants_expected_", thr, ".txt", sep="")
        write.table(ukb.sbs$SNP, out_name, quote=F, row.names=F, col.names=F, sep="\t")

        #extract dosages and merge with phenotypes
        print("    Preparing dosages...")
        dos.ph <- function.extractAndPrepare_mod(ukb.sbs, maf_thr)

        #write output: variants used for the prs
        variants <- colnames(dos.ph)
        out_name <- paste("RESULTS/PRS/variants_", thr, ".txt", sep="")
        write.table(variants[1:(length(variants)-2)], out_name, quote=F, row.names=F, col.names=F, sep="\t")

        #then calculate prs -- associations in case only at the last (p<0.5)
        print("    Calculating PRS...")
        if (thr == 0.5 & do_assoc == "yes"){ prs.res <- function.PRS_mod(dos.ph, ukb.sbs, assoc="yes") } else { prs.res <- function.PRS_mod(dos.ph, ukb.sbs, assoc="no") }
        prs <- prs.res[[1]]
        assoc <- prs.res[[2]]
        assoc.list[[c]] <- assoc

        #logit model for association
        prs$PHENO <- prs$PHENO - 1

        #now without APOE
        #first make apoe prs
        apoe.geno <- dos.ph[, c("rs7412_C", "rs429358_T", "IID", "PHENO")]
        apoe.res <- function.PRS_mod(apoe.geno, ukb.sbs, assoc="no")
        apoe.prs <- apoe.res[[1]]
        #calculate prs as difference between previous and apoe
        prs.noAPOE <- prs
        prs.noAPOE$PRS <- prs.noAPOE$PRS - apoe.prs$PRS

        #scale PRSs
        prs$PRS <- scale(x = prs$PRS)
        prs.noAPOE$PRS <- scale(x = prs.noAPOE$PRS)

        #add covariates
        covar <- read.table("GENO_DATA/20200123_covariates_kept_Age_gwas.txt", h=T, stringsAsFactors=F, sep="\t")
        prs.tmp <- merge(prs, covar, by="IID")
        prs.tmp.noAPOE <- merge(prs.noAPOE, covar, by="IID")
        
        #run model with apoe
        m <- glm(PHENO ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data=prs.tmp, family = "binomial")
        #print(summary(m))

        #run model without apoe
        m.noAPOE <- glm(PHENO ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data=prs.tmp.noAPOE, family = "binomial")
        #print(summary(m.noAPOE))

        #roc curve with APOE
        m2 <- glm(PHENO ~ PRS, data=prs.tmp, family = "binomial")
        prob <- predict(m2, type=c("response"))
        prs.tmp$prob <- prob
        g <- roc(PHENO ~ prob, data=prs.tmp)
        roc.list[[c]] <- g
        
        #roc curve without APOE
        m2.noAPOE <- glm(PHENO ~ PRS, data=prs.tmp.noAPOE, family = "binomial")
        prob <- predict(m2.noAPOE, type=c("response"))
        prs.tmp.noAPOE$prob <- prob
        g.noAPOE <- roc(PHENO ~ prob, data=prs.tmp.noAPOE)
        roc.list.noAPOE[[c]] <- g.noAPOE

        #add to general output -- with and without APOE
        model.stats[c, 1:4] <- c(length(variants), exp(summary(m)$coefficients[2, 1]), summary(m)$coefficients[2, 2], summary(m)$coefficients[2, 4])
        model.stats[c, 5:7] <- c(exp(summary(m.noAPOE)$coefficients[2, 1]), summary(m.noAPOE)$coefficients[2, 2], summary(m.noAPOE)$coefficients[2, 4])
        c <- c + 1

        #write output: actual scaled prs including APOE
        out.name <- paste("RESULTS/PRS/prs_", thr, ".txt", sep="")
        write.table(prs.tmp, out.name, quote=F, row.names=F, sep="\t")
        #write outputs: the actual scaled prss without APOE and the variants used for the prs
        out.name.noAPOE <- paste("RESULTS/PRS/prs_", thr, "_noAPOE.txt", sep="")
        write.table(prs.tmp.noAPOE, out.name.noAPOE, quote=F, row.names=F, sep="\t")

        print("    Moving to next threshold...")
}

#write stats of all models
write.table(model.stats, "RESULTS/PRS/all_models_stats_maf1pc.txt", row.names=T, sep="\t", quote=F)
write.table(assoc, "RESULTS/PRS/single_variant_all_assoc.txt", row.names=T, sep="\t", quote=F)

all.roc <- list(roc.list, roc.list.noAPOE)
save(all.roc, file="RESULTS/PRS/ROC_curve_objects.RData")

##################################

