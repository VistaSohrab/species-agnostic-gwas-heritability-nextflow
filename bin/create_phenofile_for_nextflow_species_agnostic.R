#!/usr/bin/Rscript

# Author: Vista Sohrab
# Date: 11/19/2025
# Input files: (1) A TSV file containing individual ID (IID) and phenotype columns of interest for heritability or GWAS analysis (ensure that all columns are numeric)
#              (2) A plink .fam file containing IIDs of genotyped individuals
# Goal: The goal is to create a phenotype file for GWAS and heritability analysis for each column in TSV file (except for IID)
# Output: Phenotype files generated for each column in TSV file (N phenotype files where N is the number of columns in TSV file excluding IID, therefore N-1 files will be generated)
 

# load in library
library(tidyverse)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--dataframe", 
                    help = "Tab-separated dataset of phenotypes (TSV file) for GWAS and heritability analysis")

parser$add_argument("--plinksetFam",
                    help = "Plink .fam file to identify sequenced dogs (genetic data)")

# Parse the arguments passed to the script
args <- parser$parse_args()

pheno_datafile <- args$dataframe
geneticfile <- args$plinksetFam

phenodata <- read_tsv(pheno_datafile)
 
# directly use fam file to get genotyped dog IDs
# second column of fam file contains the genotyped individual IDs (IIDs) while the first column contains FIDs set in the plinkset
sequenced_fam <- read.table(geneticfile, header = FALSE, col.names = c("FID", "IID", "father_id", "mother_id", "sex", "pheno"))


# create a dataframe of phenotypes for genotyped individuals only 
pheno_data_genotyped <- phenodata %>% filter(IID %in% sequenced_fam$IID)

extract_pheno <- function(pheno_variable){
  
  # establish first 2 columns of phenotype file (FID, and IID (dog ID))
  pheno_outfile <- sequenced_fam %>% select(FID, IID)
  
  # of the genotyped dogs with phenotypes available, choose dog id columns and phenotype column of interest
  pheno_data_genotyped_subset <- pheno_data_genotyped %>% select(IID, all_of(pheno_variable))
  
  pheno_outfile <- pheno_outfile %>% left_join(pheno_data_genotyped_subset, by = "IID")
  
  # identify number of dogs that have both phenotypes and genotypes available for genomic analysis
  phenotyped_genotyped_number <- sum(!is.na(pheno_outfile[[pheno_variable]]))
  
  # phenotype file output name
  pheno_outname <- paste0(pheno_variable, "_N-", phenotyped_genotyped_number, ".tsv")
  
  
  write.table(pheno_outfile, pheno_outname, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
}

for (pheno in colnames(phenodata)) {
  if (pheno != "IID"){
    extract_pheno(pheno)
  }
}
  
