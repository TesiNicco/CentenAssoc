# SCRIPT TO MATCH VARIANTS IN 1000G, GWAS AND OUR GENOTYPES.
# THEN ALSO PRODUCES UPDATED PLINK2 FILES THAT MATCH VARIANT IDS
# FINALLY DOES THE CLUMPING.

# LIBRARIES
library(data.table)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

# FUNCTIONS
## function to do the checkings
checkMe <- function(i, what){
    print(paste("## Start with chromosome", i))
    
    #check what needs to be done
    if (what == "match_clump"){
        print(" Reading input files..")
        # read gwas
        fname <- paste("INPUTS_OTHER/GWAS_DATA_PER_CHR/chr", i, "_gwas_data.txt", sep="")
        gwas <- fread(fname, h=T)

        # read 1000G data
        fname <- paste("GENO_DATA/LD_CLUMPING/chr", i, "_1kg.bim", sep="")
        kg <- fread(fname, h=F)

        # read my data
        fname <- paste("GENO_DATA/PLINK2/chr", i, "_plink2.pvar", sep="")
        my <- fread(fname, h=T)

        # first match is kg and gwas
        print(" Matching variants..")
        match.kg.gwas <- gwas[which(gwas$rsid %in% kg$V2),]

        # then match the matched dset with my data
        match.all <- match.kg.gwas[which(match.kg.gwas$rsid %in% my$ID),]
    
        # identify mismatches: for these try to match by chr_pos
        mism1 <- match.kg.gwas[which(!(match.kg.gwas$rsid %in% my$ID)),]
        mism2 <- my[which(!(my$ID %in% match.all$rsid)),]
        colnames(mism2) <- c("CHR", "POS", "ID", "REF", "ALT", "FILTER", "INFO")
        mism2$LOC <- paste(mism2$CHR, mism2$POS, sep="_")
        sec.match <- mism1[which(mism1$snpid %in% mism2$LOC),]
        m <- merge(sec.match, mism2, by.x="snpid", by.y="LOC")

        # merge all matches
        match.all <- rbind(match.all, sec.match)

        #exclude duplicates
        d <- m[duplicated(m$ID),]
        m <- m[which(!(m$ID %in% d$ID)), ]
        match.all <- match.all[which(!(match.all$snpid %in% d$snpid)), ]

        print(" Saving results..")
        # write variants to do clumping on
        write.table(match.all$rsid, paste("GENO_DATA/LD_CLUMPING/chr", i, "_toKeep.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
        write.table(match.all, paste("INPUTS_OTHER/GWAS_DATA_PER_CHR/chr", i, "_gwas_data.txt", sep=""), quote=F, row.names=F, col.names=T, sep="\t")
        write.table(m[, c("ID", "rsid")], paste("GENO_DATA/PLINK2/chr", i, "_convert.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")

        # run plink command to update names in my plink files
        print(" Running PLINK2 for updating names..")
        cmd <- paste("plink2 --pfile GENO_DATA/PLINK2/chr", i, "_plink2 --update-name GENO_DATA/PLINK2/chr", i, "_convert.txt --make-pgen --out GENO_DATA/PLINK2/chr", i, "_plink2_upd", sep="")
        system(cmd)

        # run plink to do the clumping
        print(" Running PLINK for clumping..")
        cmd_clump <- paste("plink --bfile GENO_DATA/LD_CLUMPING/chr", i, "_1kg --extract GENO_DATA/LD_CLUMPING/chr", i, "_toKeep.txt --maf 0.01 --clump INPUTS_OTHER/GWAS_DATA_PER_CHR/chr", i, "_gwas_data.txt --clump-field p --clump-kb 750 --clump-p1 1 --clump-r2 0.001 --clump-snp-field rsid --out RESULTS/CLUMPING/chr", i, sep="")
        system(cmd_clump)

    } else if (what == "clump_only"){
        # run plink to do the clumping
        print(" Running PLINK for clumping..")
        cmd_clump <- paste("plink --bfile GENO_DATA/LD_CLUMPING/chr", i, "_1kg --extract GENO_DATA/LD_CLUMPING/chr", i, "_toKeep.txt --maf 0.01 --clump INPUTS_OTHER/GWAS_DATA_PER_CHR/chr", i, "_gwas_data.txt --clump-field p --clump-kb 750 --clump-p1 1 --clump-r2 0.001 --clump-snp-field rsid --out RESULTS/CLUMPING/chr", i, sep="")
        system(cmd_clump)
    }
    print(paste("## Done with chromosome", i))
    return(i)
}

## function to check that everything went fine with all chromosomes
wasItRight <- function(){
    # first check plink2 files -- if everything was ok, remove previous files
    cmd = "ls GENO_DATA/PLINK2/chr*upd.pvar"
    plink2_files <- system(cmd, inter=T)
    if (length(plink2_files) != 22){ 
        print("There was a problem for a chromosome in plink2 files. Check it!") 
    } else {
        print("Good, everything went right! ")
        cmd_clean <- "for chr in {1..22}; do rm GENO_DATA/PLINK2/chr${chr}_plink2.*; done"
        system(cmd_clean)
    }

    # then check clumping files
    cmd = "ls RESULTS/CLUMPING/*clumped"
    clump_files <- system(cmd, inter=T)
    if (length(clump_files) != 22){ print("There was a problem for a chromosome in clumped files. Check it!") }

    return("Checking done!")
}

# MAIN
what <- args[1]
#help message in case
if (what == "--help"){ print("Choices are: match_clump (does matching and clumping) and clump_only (only clump)") }

res <- lapply(1:22, checkMe, what=what)
print(wasItRight)

