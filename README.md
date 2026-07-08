## Species agnostic genome-wide association study and heritability workflow 
Nextflow workflow to conduct genome-wide association studies in a given species of interest using GCTA's mixed linear model with leave-one-chromosome-out approach (GCTA MLMA LOCO) and heritability analyses using GCTA's restricted maximum likelihood


### Input files 

The input files include the following:

  * A phenotype file containing the column `IID` and columns for each phenotype of interest to be analyzed (tab-delimited format with headers: IID"\t"phenoA"\t"phenoB)

```
IID  F1      F2      F3      F4      F5      F6      F7      F8
28      5.114824997780695       4.11126844119854        -1.367657516345442      3.2407159840647464      15.110542318327791      3.252383263718634       -8.958823306743273      -0.8692999350381554
46      12.041305639988916      10.18177012301823       -1.415599633101153      0.5474706316540905      17.082930711862293      9.72219562232409        -11.800388098421916     -5.28880991642094
95      9.786001103329315       16.189784949097383      -1.910598665812524      -11.137085834669682     23.089078417388883      -1.8888784705935753     -10.38988615567527      -5.28880991642094
98      1.884329900709555       NA      NA      NA      NA      NA      NA      NA
```

  * Genetic plinkset in PLINK1 (.bim/.bed/.fam) format (filtering of the plinkset occurs within the workflow)
  
  * Reference genome annotation file in the following format (tab-delimited file containing chr, start position, end position, gene name)

 ```
1       13565063        13577444        SERPINB13
1       13587694        13611747        SERPINB12
1       13639168        13716725        SERPINB5
1       13671047        13674721        LOC119870705
1       13716837        13751163        VPS4B
1       13766594        13804675        KDSR
1       13811586        13978482        BCL2
```
#### Optional input files 

  * A quantitative covariate file containing the columns `IID`, `FID`, and columns for each quantitative covariate of interest (such as age and weight) but without header in tab-delimited format
    
```
0       28      1.68407 3
0       46      1.5219800000000001      3
0       95      7.53571 3000000001      3
0       98      13.8187 2
```

  * A discrete covariate file containing the columns `IID`, `FID`, and columns for each discrete covariate of interest (such as sex) but without header in tab-delimited format

```
0       28      female
0       46      male
0       95      female
0       98      female  
```

### Parameters 
Covariate files are optional. To add a discrete covariate, provide full path to discrete covariate file of interest with --covar and path to quantitative covariate file with --qcovar.
Chr X is excluded from analyses by default, but if you would like to include it in analyses, then set --include_chrX = "true" at job submission.
Default values are listed that should be changed upon submission of workflow especially autosome number depending on the species (default is set for dogs).
Path to input files are required for --geneticset, --annotation, and --pheno_dateset.

```
params.geneticset = "/scratch/vsohrab/darwins_dogs/darwins_dogs_genetic_set_2024/DarwinsDogs_2024_N-3277_canfam4_gp-0.70_biallelic.{bed,bim,fam}"

params.maf = 0.01 (minor allele frequency threshold - default value is 0.01)

params.geno = 0.05 (genotyping threshold - default value is 0.05)

params.hwe = 1e-20 (Hardy-Weinberg equilibrium p-value - default value is 1e-20)

params.autosome_num = 38 (number of autosomes - default value is set to 38 for dogs)

params.ld_window_kb = 250 (specifying the length of the region for segment-based LD score calucation of GCTA - default value is 250 kb) 

params.clump_pval = 1e-6 (significance threshold for index SNPs - default value is 1e-6)

params.clump_kb = 250 (physical distance threshold for clumping
 - default value is 250 kb)

params.clump_rsquared = 0.2 (LD threshold for clumping - default is 0.2)

params.annotation = "/scratch/vsohrab/reference/UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat.genes.validchr.bed"

params.pheno_dataset = /scratch/vsohrab/dog_gwas/DarwinsArk_8Factors_N-3277.tsv

params.covar = null

params.qcovar = null 

params.include_chrX = false

```

#### Note
Each workflow run is unique to a specific discrete covariate and quantitative covariate set pair (for example, if I have 10 phenotypes to analyze where quantitative covariates are age and weight and discrete covariate is sex; then, I can include those 10 phenotypes in a single phenotype input file with IID, pheno1, pheno2, pheno3,...,pheno10; however, if I have another set of phenotypes with a different set of covariates or would like to exclude covariates altogether, then I would need to run those analyses separately. If quantitative covariate is age and discrete covariate is sex, then I will create a separate phenotype file for that workflow run with IID, phenoA, phenoB). If running the same phenotype with and/or without a particular covariate, embed the covar and qcovar name in the phenotype input file column name, because the column names in the input phenotype file is used to create the output files (and if those are the same names, then the will be overwritten with the most recent workflow run when running the same phenotype with different set of covariates in the same working directory)


### Output files
* Filtered genetic plinkset used in analyses can be found in *filtered_plinkset*
* The full genetic relatedness matrix calculated for analyses can be found in *grm*
* GWAS output file(s) are reported in *association* folder
* Heritability output file(s) are reported in *reml* folder
* Phenotype file(s) generated by workflow can be found in *phenofiles* folder
* Clump output file(s) are reported in *clump_pval-$pval_kb-$kb_r2-$r2* folder. For example if all default parameters are used, clump output files can be found in *clump_pval-0.000001_kb-250_r2-0.2*
* Manhattan plots and qqplots are generated in *plots* folder
* For viewing an example log file from GWAS result, please check the *logs* folder (ie to view sample size of GWAS analysis for a particular phenotype)


### Example code to run genetic analysis: 

with discrete covariate file only: 
```
/PATH/TO/nextflow run /PATH/TO/species_agnostic_gwas_heritability_multiple_phenotypes.nf \
-c /PATH/TO/nextflow_slurm_unity_dog_gwas_heritability.config \
-resume -w /PATH/TO/work \
--geneticset "/PATH/TO/genetic_set_prefix.{bed,bim,fam}" \
--pheno_dataset /PATH/TO/phenotype_input_file.tsv \
--covar /PATH/TO/discrete_covar_input_file.tsv \
--annotation /PATH/TO/annotation.bed \
--outdir /PATH/TO/workdir \
--maf ${maf} --geno ${geno} --hwe ${hwe}
```

with quantitative covariate file only and including chrX in association analyses: 
```
/PATH/TO/nextflow run /PATH/TO/species_agnostic_gwas_heritability_multiple_phenotypes.nf \
-c /PATH/TO/nextflow_slurm_unity_dog_gwas_heritability.config \
-resume -w /PATH/TO/work \
--geneticset "/PATH/TO/genetic_set_prefix.{bed,bim,fam}" \
--pheno_dataset /PATH/TO/phenotype_input_file.tsv \
--qcovar /PATH/TO/qcovar_input_file.tsv \
--annotation /PATH/TO/annotation.bed \
--outdir /PATH/TO/workdir \
--maf ${maf} --geno ${geno} --hwe ${hwe} \
--include_chrX "true"
```


This is an updated version of https://github.com/VistaSohrab/dog-gwas-heritability-nextflow. With thanks to Kathleen Morrill (for Manhattan and qqplot R scripts) and Rob Bierman for useful discussion in developing the first version of workflow. 


