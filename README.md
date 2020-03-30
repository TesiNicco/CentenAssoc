# CentenAssoc
Pipeline to perform replication of genetic variants associated with longevity in previous studies by looking at the reported variants, their haplotype, the genes these variants are associated with, and polygenic risk scores.

## Set-up working environment
First and (perhaps) most important, we need to set-up the environment with the necessary folders and files. To guide the reader through the set-up, we will do folder-by-folder:

### GENO_DATA folder
1. Create a folder for the genotype data (GENO_DATA): please put the genotype files (PLINK1.9 .bed/.bim/.fam files) in this folder, together with the relative phenotype, sample and covariates file. The format of the latter three files should be as reported in PLINK website (https://www.cog-genomics.org/plink/1.9/)
2. In the GENO_DATA folder, creates a new folder for the PLINK2 (PLINK2) files: in principle, this is not strictly necessary as the PRS can also be calculated with the PLINK .bim/.bed/.fam files. However, we noticed that missing values were present when extracting dosages (via the --export A option in PLINK1.9), while they were not present when using PLINK2 .pgen/.pvar/.psam files
3. Create a folder where to put the scripts (BIN): scripts will be downloaded from Github (this page) with the close command
4. Finally, also create folders for the results (RESULTS) and for additional inputs required (INPUTS_OTHER)
At the end of this step, you should have your envirnment set. It should include the pipeline file (pipeline.sh) and 4 folders (BIN, RESULTS, GENO_DATA and INPUTS_OTHER).

## Overview of data format files
Before you start, please make sure that all the files and tools required are in place. The pipeline requires R installed and working in your system. Additional tools are based on Unix sheel (awk, mkdir, paste, cut). Therefore, our pipeline is tested and workig on Unix systems only.  

In addition, the pipeline in the current stage assumes that the name of the files are as it follows:
1. genotype files: chrXX_dose.bim/.bed/.fam (in GENO_DATA/ folder)
2. phenotype file: 20200317_phenotypes_subsetAge.txt
3. samples file: 20200317_samples_subsetAge.txt
4. covariate file: 20200123_covariates_kept_Age_gwas.txt
The format of these files should be as indicated in PLINK documentation (https://www.cog-genomics.org/plink/1.9/)
5. 
