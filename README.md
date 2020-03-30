# CentenAssoc
Pipeline to perform replication of genetic variants associated with longevity in previous studies by looking at the reported variants, their haplotype, the genes these variants are associated with, and polygenic risk scores.

## Set-up working environment
First and (perhaps) most important, we need to set-up the environment with the necessary folders and files. First thing to do is to get the git and all provided files.
```
wget https://github.com/TesiNicco/CentenAssoc.git .
```


### GENO_DATA folder
```
mkdir GENO_DATA
mkdir GENO_DATA/PLINK2
mkdir GENO_DATA/LD_CLUMPING
```
Place in the `GENO_DATA/` folder:
- Genotype files (PLINK1.9 .bed/.bim/.fam files), together with the relative phenotype, sample and covariates file. The format of the latter three files should be as reported in PLINK website (https://www.cog-genomics.org/plink/1.9/)
- For the survival analysis, data regarding the survival are required: these include the APOE status of all individuals and the phenotypes of the individuals with follow-up data for the survival analysis. 

Place in the `GENO_DATA/PLINK2/` folder:
- Genotype files (PLINK2 .pgen/.pvar/.psam files): in principle, this is not strictly necessary as the PRS can also be calculated with the PLINK .bim/.bed/.fam files. However, we noticed that missing values were present when extracting dosages (via the --export A option in PLINK1.9)

Place in the `GENO_DATA/LD_CLUMPING/` folder:
- Genotype files of the 1000Genome project (Phase 3), european samples only. You can download them from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/. Selecting only European sample will be important in order to match population. You can look up the samples information and download them from https://www.internationalgenome.org/data-portal/sample. Once data is downloaded, you should convert the VCF files into PLINK. We suggest to create PLINK files with european individuals only, and restrict to variants with minor allele frequency of 1% or higher. This will greatly reduce the size of the genotype files. Strictly speaking, using genotypes from an independent cohort is highly suggestable but not mandatory. In principle, these genotypes are used for the calculation of LD and LD-clumping. You can (although we do not encourage) use your genotypes here.

### BIN folder
```
mkdir BIN/
cd BIN/
wget https://github.com/TesiNicco/CentenAssoc.git .
```
In this folder all the scripts will be placed.

### INPUTS_OTHER folder
```
mkdir INPUTS_OTHER/
mkdir INPUTS_OTHER/CADD/
mkdir INPUTS_OTHER/GWAS_DATA_PER_CHR/
mkdir INPUTS_OTHER/RECOMB_RATES/
```
This folder includes all annotation files that will be used in the analyses. You need to download these files and place them in the correct folders. Please place in `INPUTS_OTHER/`
- Ensemble_to_GeneName.txt: file to map ENSEMBLE IDs to gene name (provided in the git)
- GWAS_catalog_20200310.txt: file downloadable from GWAS catalog (https://www.ebi.ac.uk/gwas/api/search/downloads/full)
- GWAS_data_Timmers_cleaned.txt: summary statistics of the GWAS by-proxy on longevity (Timmers et al., 2019, eLife) available at https://datashare.is.ed.ac.uk/handle/10283/3209
- hg19_geneListandPos: gene names and positions from RefSeq (provided in the git)
- NCBI37.3.gene.loc: gene positions to use with MAGMA (provided in the git)
- reported_variants.txt: chromosome, positions and IDs of the reported variants to test (provided in the git)
- Whole_Blood.v7.signif_variant_gene_pairs.txt.gz: GTEx dataset of variant-gene associations from the blood (downloadable at https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz)

Please place in `INPUTS_OTHER/GWAS_DATA_PER_CHR/`:
- the summary statistics file (GWAS_data_Timmers_cleaned.txt) divided by chromosomes

Please place in `INPUTS_OTHER/CADD/`:
- CADD annotation files: download the annotation file of the 1000Genome project including annotation (1000 phase 3 Genome variants (SNVs and InDels) incl. all annotations, 6.1GB) from https://cadd.gs.washington.edu/download. Once downloaded, to store items more efficiently and boost computations, you should divide the main file in 22 chromosome files. You can do this with bash commands (for `example` for loop and cut command). The resulting 22 files should be placed in the folder and named (chrXX_cadd.txt). For your information, only columns `#Chrom Pos Ref Alt Type AnnoType Consequence GeneName PHRED` will be used in the script, if you want to restrict to these column only when dividing whole file in chromosome-files.

Please place in `INPUTS_OTHER/RECOMB_RATES/`:
- recombination rates for each chromosome (provided in the git)
