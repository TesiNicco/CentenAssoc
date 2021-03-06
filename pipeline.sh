##############################################
# FIRST, SET UP ENVIRONMENT FOR ALL ANALYSES #
# GENERAL DEFINITION OF CASES: CENTENARIANS  #
# GENERAL DEFINITION OF CONTRONLS:	     #
# COMBINATION OF: LASA, SCD, BRAIN BANK,     #
# TWIN STUDY AND 100-PLUS CONTROLS	     #
#                                            #
#                                            #
# FOR SIMPLICITY, ALL STEPS REQUIRING RAW    #
# GENOTYPE FILES AND LARGE DATA ARE MARKED   #
# (***) AS THESE FILES WILL NOT BE PUBLICLY  #
# AVAILABLE BUT THEY CAN BE REQUESTED.       #
##############################################

##########################################################################
# CREATE FOLDERS FOR RESULTS                                             #
##########################################################################
mkdir RESULTS

##########################################################################
# DONE WITH SETTING UP THE ENVIRONMENT                                   #
##########################################################################

##########################################################################
# START WITH THE ANALYSES                                                #
# FIRST, SINGLE-VARIANT ANALYSIS                                         #
##########################################################################
## step 1 -- create directory for results
mkdir RESULTS/SINGLE_VARIANT/
## step 2 -- given the variants to look at, check for id for all of them
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
while read -r chrom pos id; do grep -w ${pos} GENO_DATA/chr${chrom}_dose.bim | cut -f2 >> RESULTS/SINGLE_VARIANT/reported_variants_ids.txt ; done <INPUTS_OTHER/reported_variants.txt
while read -r chrom pos id; do grep -w ${pos} GENO_DATA/chr${chrom}_dose.bim | sed 's/ /\t/g' >> RESULTS/SINGLE_VARIANT/reported_variants_info.txt ; done <INPUTS_OTHER/reported_variants.txt
## step 3 -- do associations, blindly on all chromosomes and the selected variants
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
for chr in {1..22}; do plink2 --bfile GENO_DATA/chr${chr}_dose --covar GENO_DATA/20200123_covariates_kept_Age_gwas.txt --covar-name PC1,PC2,PC3,PC4,PC5 --extract RESULTS/SINGLE_VARIANT/reported_variants_ids.txt --mach-r2-filter 0.8 2.0 --glm hide-covar cols=+beta,+a1freq,+a1freqcc,+machr2 --out RESULTS/SINGLE_VARIANT/chr${chr}_tmp_assoc; done
## step 4 -- put stuff together and clean
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
mv RESULTS/SINGLE_VARIANT/chr1_tmp_assoc.PHENO1.glm.logistic RESULTS/SINGLE_VARIANT/chrAll_reported_variants_assoc.txt
for chr in {2..22}; do grep -v "BETA" RESULTS/SINGLE_VARIANT/chr${chr}_tmp*logistic >> RESULTS/SINGLE_VARIANT/chrAll_reported_variants_assoc.txt ; done
rm RESULTS/SINGLE_VARIANT/*logistic
rm RESULTS/SINGLE_VARIANT/*log
Rscript BIN/checkReported.R

##########################################################################
# GENE-BASED ANALYSIS -- this will run MAGMA tool, so make sure that is  #
# correctly installed and working in your system.			 #
# In addition, to test genes, we will need to map variants to genes.     #  
##########################################################################
## step 1 -- create directory for results and for annotation
mkdir RESULTS/GENE_BASED
mkdir RESULTS/VARIANT_GENE_MAPPING
## step 2 -- perform variant-gene mapping, that will generate annotation file and extract MAGMA coordinates
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
## The whole CADD is required and it is not provided here. The same analysis can be run through snpXplorer (http://snpxplorer.eu.ngrok.io) using the reported variants as input for AnnotateMe. Alternatively download the whole CADD and use the following script
## Rscript BIN/annotate_reported.R
## step 3 -- need to modify covariates for MAGMA 
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
awk 'NR>1 {print 0, $1, $2, $3, $4, $5, $6}' GENO_DATA/20200123_covariates_kept_Age_gwas.txt | sed 's/ /\t/g' > RESULTS/GENE_BASED/covariates_chc_ctr.tab
## step 4 -- need to calculate inflation and annotations, so create first folders to make order
mkdir RESULTS/GENE_BASED/INFLATION
mkdir RESULTS/GENE_BASED/ANNOTATIONS
## step 5 -- then randomly take 5000 genes to calculate inflation
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
sort -r INPUTS_OTHER/NCBI37.3.gene.loc | head -5000 > RESULTS/GENE_BASED/INFLATION/random_5k_genes_inflation.txt
## step 6 -- perform MAGMA annotation
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
for chr in {1..22}; do magma --annotate window=2 --snp-loc GENO_DATA/chr${chr}_dose.bim --gene-loc RESULTS/GENE_BASED/INFLATION/random_5k_genes_inflation.txt --out RESULTS/GENE_BASED/INFLATION/chr${chr}_5k; done
## step 7 -- run MAGMA test
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
for chr in {1..22}; do magma --bfile GENO_DATA/chr${chr}_dose --genes-only --gene-annot RESULTS/GENE_BASED/INFLATION/chr${chr}_5k.genes.annot --covar file=RESULTS/GENE_BASED/covariates_chc_ctr.tab --gene-model multi=snp-wise --out RESULTS/GENE_BASED/INFLATION/chr${chr}_5k_assoc; done
## step 8 -- merge inflation results
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
for chr in {2..22}; do grep -v GENE RESULTS/GENE_BASED/INFLATION/chr${chr}_5k_assoc.genes.out >> GENE RESULTS/GENE_BASED/INFLATION/chr1_5k_assoc.genes.out; rm GENE RESULTS/GENE_BASED/INFLATION/chr${chr}_5k_assoc.genes.out ; done
mv RESULTS/GENE_BASED/INFLATION/chr1_5k_assoc.genes.out RESULTS/GENE_BASED/INFLATION/chrAll_5k_assoc.genes.out
## step 9 -- look at inflation factor
Rscript BIN/calculateInflation.R RESULTS/GENE_BASED/INFLATION/chrAll_5k_assoc.genes.out
## step 10 -- if inflation if fine, go to the analysis of actual gene-set -- first, annotation step
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
for chr in {1..22}; do magma --annotate window=2 --snp-loc GENO_DATA/chr${chr}_dose.bim --gene-loc RESULTS/GENE_BASED/gene_coordinates.txt --out RESULTS/GENE_BASED/ANNOTATIONS/chr${chr}; done
## step 11 -- do actual MAGMA test
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
for chr in {1..22}; do magma --bfile GENO_DATA/chr${chr}_dose --genes-only --gene-model multi=snp-wise --gene-annot RESULTS/GENE_BASED/ANNOTATIONS/chr${chr}.genes.annot --covar file=RESULTS/GENE_BASED/covariates_chc_ctr.tab --out RESULTS/GENE_BASED/chr${chr}_multiSNP_2kb_annot; done
## step 12 -- put results together
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
for chr in {2..22}; do grep -v START RESULTS/GENE_BASED/chr${chr}_multiSNP_2kb_annot.genes.out >> RESULTS/GENE_BASED/chr1_multiSNP_2kb_annot.genes.out; rm RESULTS/GENE_BASED/chr${chr}_multiSNP_2kb_annot.genes.out; done
mv RESULTS/GENE_BASED/chr1_multiSNP_2kb_annot.genes.out RESULTS/GENE_BASED/chrAll_multiSNP_2kb_annot.genes.out
## step 13 -- correct pvalues and analyze results
Rscript BIN/analyzeMagma.R

##########################################################################
# PRS ANALYSIS -- this will first run the clumping with PLINK and then   #
# will run a script to perform the PRS at different p-value thresholds.  #
# The script will also produce some plots and the survival analysis.     #
##########################################################################
## step 1 -- create directories for clumping, PRS and genotypes for LD
mkdir RESULTS/CLUMPING
mkdir RESULTS/PRS
mkdir RESULTS/PRS_ANNOTATION
## step 2 -- match variants between 1000G, GWAS and our genotypes and do the clumping
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
Rscript BIN/matchAndClump.R [remember to add an option, either: match_clump OR clump_only]
## step 3 -- put clumped variants together
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
echo "CHR SNP BP P TOTAL" | sed 's/ /\t/g' > RESULTS/CLUMPING/clumping_chrAll_variantsKept.txt
for chr in {1..22}; do awk '{print $1, $3, $4, $5, $6}' RESULTS/CLUMPING/chr${chr}.clumped | sed 's/ /\t/g' | grep -v CHR >> RESULTS/CLUMPING/clumping_chrAll_variantsKept.txt; done
## step 4 -- calculate PRS on the reported variants and at different thresholds
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
Rscript BIN/build_PRS.R [yes/no for single variant association]
## step 5 -- make some plots, do the survival analysis with the best model
## *** OUTPUT FILES ARE PLACED IN THE RESPECTIVE DIRECTORIES
Rscript BIN/PRS_plots_and_survival.R
## step 6 -- perform functional annotation and GWAS catolog matching
## perhaps use AnnotateMe for this
## The whole CADD is required and it is not provided here. The same analysis can be run through snpXplorer (http://snpxplorer.eu.ngrok.io) using the PRS variants as input for AnnotateMe. Alternatively download the whole CADD and use the following script
Rscript BIN/funcAnnot_GWAScat.R
## step 7 -- gene expression
Rscript geneExpression.R

##########################################################################
# PLOTS -- figure 1 should be the regional plots of the FDR significant  #
# associations: this will comprise APOE, ABO, GBX2 and CDKN2B region.    #
# In order to make these plots, I need the single variant associations   #
# of all variants in the regions. Also remember to put the files for     #
# the recombination rates in the INPUTS_OTHER/RECOMB_RATES folder	 #
##########################################################################
mkdir RESULTS/REGIONAL_PLOTS
mkdir INPUTS_OTHER/RECOMB_RATES
Rscript BIN/RegionalPlots.R
