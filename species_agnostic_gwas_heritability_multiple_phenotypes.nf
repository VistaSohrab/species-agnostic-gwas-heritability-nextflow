#!/usr/bin/env nextflow

params.geneticset = '/seq/vgb/vsohrab/gwas_heritability_nextflow_pipeline/*{bed,bim,fam}'
params.maf = 0.01
params.geno = 0.05
params.hwe = 1e-20
params.autosome_num = 38
params.ld_window_kb = 250
params.clump_pval = 1e-6
params.clump_kb = 250
params.clump_rsquared = 0.2
params.annotation = '/seq/vgb/dap/gwas/annotate/Canis_lupis_familiaris.CanFam3.1.ensembl.gene_annotatations.withHuman.bed'
params.pheno_dataset='/seq/vgb/vsohrab/dap/gwas/pheno/DAP_2022_baseline_HLES_AFUS_CSLB_GeneticData_Dec2022_N6358/scaled_cleaned_DAP2022_baseline_survey_20231219.tsv'
params.covar = null
params.qcovar = null
params.include_chrX = false


process FILTER_PLINKSET {

    publishDir (params.outdir + '/filtered_plinkset'), mode: 'copy'

    input: 
    
    tuple (val(plinkset_prefix), path(plinkfiles))
    val maf
    val geno
    val hwe
    val autosome_num

    output:

    tuple (val("${plinkset_prefix}_maf${maf}_geno${geno}_hwe${hwe}"), path("${plinkset_prefix}_maf${maf}_geno${geno}_hwe${hwe}.{bed,bim,fam}"))

    script:

    """
    /usr/bin/plink2/plink2 \
       --bfile ${plinkset_prefix} \
       --chr-set ${autosome_num} \
       --maf ${maf} \
       --geno ${geno} \
       --hwe ${hwe} midp keep-fewhet \
       --make-bed \
       --out "${plinkset_prefix}_maf${maf}_geno${geno}_hwe${hwe}"

    """
}

process CREATE_FULL_GRM {

    publishDir (params.outdir + '/grm'), mode: 'copy' 

    input:

    tuple (val(filtered_plinkset_prefix), path(filtered_plinkfiles))
    val autosome_num

    output:

    tuple (val(filtered_plinkset_prefix), path("${filtered_plinkset_prefix}.{grm.id,grm.bin,grm.N.bin}"))

    script:

    """
    /usr/bin/gcta-1.94.1/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
       --bfile ${filtered_plinkset_prefix} \
       --make-grm \
       --autosome \
       --autosome-num ${autosome_num} \
       --out ${filtered_plinkset_prefix} \
       --thread-num 8
    """

}


process CALCULATE_CHRLEVEL_LDSCORE {

    input:

    tuple (val(filtered_plinkset_prefix), path(filtered_plinkfiles))
    val ld_window_kb
    val autosome_num
    each chr

    output:

    path "${filtered_plinkset_prefix}_kb${ld_window_kb}_chr${chr}.score.ld"

    script:

    """   
    /usr/bin/gcta-1.94.1/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
       --bfile ${filtered_plinkset_prefix} \
       --ld-score-region ${ld_window_kb} \
       --autosome-num ${autosome_num} \
       --chr ${chr} \
       --out "${filtered_plinkset_prefix}_kb${ld_window_kb}_chr${chr}" \
       --thread-num 8
    """

}

process GENERATE_PHENOFILE {
  
         
    publishDir (params.outdir + '/phenofiles'), mode: 'copy'

    input: 
    path pheno_dataset
    tuple (val(plinkset_prefix), path(plinkfiles))
    
    output:

    path("*.tsv")

    script:
    
    """
    create_phenofile_for_nextflow_species_agnostic.R --dataframe ${pheno_dataset} --plinksetFam ${plinkset_prefix}.fam

    """

}

process CREATE_LDSTRATIFIED_GRMS {

    input:
    
    path(chrlevel_ldscore)
    tuple (val(filtered_plinkset_prefix), path(filtered_plinkfiles))
    val autosome_num

    output:

    path "paths_ldstratified_mgrmlist.txt"

    script:

    """
    ld_score_quartiles.R

    for quartile in 'Q000-Q025' 'Q025-Q050' 'Q050-Q075' 'Q075-Q100'
    do
   
       /usr/bin/gcta-1.94.1/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
          --bfile ${filtered_plinkset_prefix} \
	  --make-grm \
	  --autosome-num ${autosome_num} \
	  --thread-num 8 \
	  --extract "merged.score.ld.\${quartile}.txt" \
	  --out "${filtered_plinkset_prefix}.\${quartile}"
        
       echo "\$PWD/${filtered_plinkset_prefix}.\${quartile}" >> paths_ldstratified_mgrmlist.txt
     done   

    """
}


process CALCULATE_HERITABILITY_REML_LDSTRATIFIED_CONSTRAINED {

     publishDir (params.outdir + '/reml'), mode: 'copy'

     errorStrategy 'ignore'


     input: 

     path mgrmlist
     each phenofile
     path covarfile
     path qcovarfile

     output:

     path "${phenofile.baseName}.REML.lds.hsq"

     script:   
     def covar_args = []

     // Conditionally add arguments to the list
     if (params.qcovar) {
         covar_args.add("--qcovar ${qcovarfile}")
     }

     if (params.covar) {
         covar_args.add("--covar ${covarfile}")
     }

     def optional_args = covar_args.join(' ')

     """
     
     /usr/bin/gcta-1.94.1/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
         --mgrm ${mgrmlist} \
	 --reml \
	 --pheno ${phenofile} \
	 --mpheno 1 \
	 ${optional_args} \
	 --thread-num 8 \
	 --out "${phenofile.baseName}.REML.lds"

     """
}

process CALCULATE_HERITABILITY_REML_LDSTRATIFIED_NOCONSTRAINT {

     publishDir (params.outdir + '/reml'), mode: 'copy'

     errorStrategy 'ignore'


     input:

     path mgrmlist
     each phenofile
     path covarfile
     path qcovarfile

     output:

     path "${phenofile.baseName}.REML.lds.no-constraint.hsq"

     script:

     def covar_args = []

     // Conditionally add arguments to the list
     if (params.qcovar) {
         covar_args.add("--qcovar ${qcovarfile}")
     }

     if (params.covar) {
         covar_args.add("--covar ${covarfile}")
     }

     def optional_args = covar_args.join(' ')

     """
     
     /usr/bin/gcta-1.94.1/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
       --mgrm ${mgrmlist} \
       --reml-no-constrain \
       --pheno ${phenofile} \
       --mpheno 1 \
       ${optional_args} \
       --thread-num 8 \
       --out "${phenofile.baseName}.REML.lds.no-constraint"

     """

}


process CALCULATE_HERITABILITY_REML_CONSTRAINED {

     publishDir (params.outdir + '/reml'), mode: 'copy'

     errorStrategy 'ignore'

     input:
     
     tuple (val(grm_basename), path(grmfiles))
     each phenofile
     path covarfile
     path qcovarfile

     output:
     
     path "${phenofile.baseName}.REML.no-lds.hsq"


     script:
     
     def covar_args = []

     // Conditionally add arguments to the list
     if (params.qcovar) {
         covar_args.add("--qcovar ${qcovarfile}")
     }

     if (params.covar) {
         covar_args.add("--covar ${covarfile}")
     }

     def optional_args = covar_args.join(' ')


     """
     
     /usr/bin/gcta-1.94.1/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
        --grm ${grm_basename} \
        --reml \
        --pheno ${phenofile} \
        --mpheno 1 \
        ${optional_args} \
        --out "${phenofile.baseName}.REML.no-lds"
 
    """

}


process CALCULATE_HERITABILITY_REML_NOCONSTRAINT {

          
     publishDir (params.outdir + '/reml'), mode: 'copy'
    
     errorStrategy 'ignore'

     input:
     
     tuple (val(grm_basename), path(grmfiles))
     each phenofile
     path covarfile
     path qcovarfile

     output: 

     path "${phenofile.baseName}.REML.no-lds.no-constraint.hsq"


     script:
     def covar_args = []

     // Conditionally add arguments to the list
     if (params.qcovar) {
         covar_args.add("--qcovar ${qcovarfile}")
     }

     if (params.covar) {
         covar_args.add("--covar ${covarfile}")
     }

     def optional_args = covar_args.join(' ')


     """
     /usr/bin/gcta-1.94.1/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
        --grm ${grm_basename} \
        --reml-no-constrain \
        --pheno ${phenofile} \
        --mpheno 1 \
        ${optional_args} \
        --out "${phenofile.baseName}.REML.no-lds.no-constraint"

     """
}


process GWAS_MLMA_LOCO {
     
     input:

     tuple (val(filtered_plinkset_prefix), path(filtered_plinkfiles))
     tuple (val(grm_basename), path(grmfiles))
     each chr
     val autosome_num
     each phenofile
     path covarfile
     path qcovarfile
    
     output:

     tuple(val("${phenofile.baseName}"), path("${phenofile.baseName}_${chr}.mlma"), path("${phenofile.baseName}_${chr}.log"))

     script:
    
     def covar_args = []

     // Conditionally add arguments to the list
     if (params.qcovar) {
         covar_args.add("--qcovar ${qcovarfile}")
     }

     if (params.covar) {
         covar_args.add("--covar ${covarfile}")
     }

     def optional_args = covar_args.join(' ')



     """
     if [[ ${chr} -le ${autosome_num} ]]
     then
           /usr/bin/gcta-1.94.1/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
              --bfile ${filtered_plinkset_prefix} \
              --make-grm \
              --autosome-num ${autosome_num} \
              --chr ${chr} \
              --out ${filtered_plinkset_prefix}_chr${chr} \
              --thread-num 8

           /usr/bin/gcta-1.94.1/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
              --bfile ${filtered_plinkset_prefix} \
              --grm ${grm_basename} \
              --mlma \
              --mlma-subtract-grm ${filtered_plinkset_prefix}_chr${chr} \
              --autosome-num ${autosome_num} \
              --chr ${chr} \
              --thread-num 8 \
              --pheno ${phenofile} \
              ${optional_args} \
              --out ${phenofile.baseName}_${chr}
     else 
         
           /usr/bin/plink/plink \
             --bfile ${filtered_plinkset_prefix} \
             --chr-set ${autosome_num} \
             --keep-allele-order \
             --allow-no-sex \
             --chr X \
             --make-bed \
             --out ${filtered_plinkset_prefix}_chr${chr}

       
           /usr/bin/gcta-1.94.1/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
              --bfile ${filtered_plinkset_prefix}_chr${chr} \
              --grm ${grm_basename} \
              --mlma \
              --thread-num 8 \
              --pheno ${phenofile} \
              --mpheno 1 \
              ${optional_args} \
              --out ${phenofile.baseName}_${chr}
     fi

     """
}

process MERGE_GWAS_RESULTS {

     publishDir (params.outdir + '/association'), mode: 'copy', pattern:'*.loco.mlma'
     publishDir (params.outdir + '/logs'), mode: 'copy', pattern: '*.log'

     input: 

     tuple (val(pheno), path(per_chr_gwas_results), path(per_chr_gwas_logfiles))

     output:

     tuple(val("${pheno}"), path("${pheno}.loco.mlma"), path("${per_chr_gwas_logfiles[0]}"))

     script:

     """
     touch "${pheno}.loco.mlma"

     # extract header from first mlma loco output file (typically header from chr1 file)
     head -n1 ${per_chr_gwas_results[0]} > "${pheno}.loco.mlma"

     for val in ${per_chr_gwas_results}
     do
         # merge chr level mlma loco output files together into a single .loco.mlma file
         tail -n+2 \${val} >> "${pheno}.loco.mlma"
     done

     """
}


process CLUMP {

    
     publishDir (params.outdir + '/clump_pval-' + params.clump_pval + '_kb-' + params.clump_kb + '_r2-' + params.clump_rsquared), mode: 'copy', pattern: '*.clumped.ranges'

     errorStrategy 'ignore'

     input:

     tuple (val(outname), path(phenofile), path(gwas_mlma_loco_outfile), path(gwas_mlma_loco_logfile), val(filtered_plinkset_prefix), path(filtered_plinkfiles))
     val clump_pval
     val clump_kb
     val clump_r2
     val autosome_num
     path annotation

 
     output:
     
     tuple (path("${outname}.clumped"), path ("${outname}.clumped.ranges"))

     script:

     """
     /usr/bin/plink/plink \
       --bfile ${filtered_plinkset_prefix} \
       --clump ${gwas_mlma_loco_outfile} \
       --allow-no-sex \
       --chr-set ${autosome_num} \
       --pheno ${phenofile} \
       --clump-field p \
       --clump-p1 ${clump_pval} \
       --clump-kb ${clump_kb} \
       --clump-r2 ${clump_r2} \
       --clump-range ${annotation} \
       --out ${outname}

     """
}

process PLOT_MANHATTAN_QQPLOT {

     publishDir (params.outdir + '/plots'), mode: 'copy'
     
     errorStrategy 'ignore'

     input:
     
     tuple (val(outname), path(clumped_outfile), path(gwas_mlma_loco_outfile), path(gwas_mlma_loco_logfile))

     output:
     
     path "${outname}.manhattan.png"
     path "${outname}.qqplot.png"

     script: 

     """
     manhattan_plot.R --mlma ${gwas_mlma_loco_outfile} --clump ${clumped_outfile} --title ${outname} --outname ${outname}

     """

}

workflow {

     
     //create a channel of plink files from raw genetic set (transfer .bim, .bed, .fam to FILTER_PLINKSET in workflow)
     raw_plinkset_ch = Channel.fromFilePairs(params.geneticset, size:-1)

     //create a channel of autosome numbers for LD score calculation
     autosome_ch = Channel.of(1..params.autosome_num)

     //create a channel of autosome numbers up to autosome_num + 1 for GWAS (where autosome_num + 1 signifies X chr) if indicated by include_chrX
     gwas_chr_ch = params.include_chrX ? Channel.of(1..(params.autosome_num+1)) : Channel.of(1..params.autosome_num)
     
    //define covariate and quantitative covariate optional input parameters
     def covar = params.covar ? file(params.covar) : []
     def qcovar = params.qcovar ? file(params.qcovar) : []
     
    // filter raw genetic plinkset based on minor allele frequency, genotyping rate, and hardy-weinberg equilibrium p-value threshold
     FILTER_PLINKSET(raw_plinkset_ch, params.maf, params.geno, params.hwe, params.autosome_num)


     // create complete grm using filtered genetic plinkset and provide the number of autosomes (ie for dogs: 38)
     CREATE_FULL_GRM(FILTER_PLINKSET.out, params.autosome_num)
 
     // calculate LD scores for SNPs at the chr level - inputs are filtered plinkset and autosome number
     CALCULATE_CHRLEVEL_LDSCORE(FILTER_PLINKSET.out, params.ld_window_kb, params.autosome_num, autosome_ch)

     // gather all resulting chr level LD score files
     all_ldscore_results = CALCULATE_CHRLEVEL_LDSCORE.out.collect()

     // generate phenotype files from phenotype dataset TSV file     

     GENERATE_PHENOFILE(params.pheno_dataset, raw_plinkset_ch)
      
    
     // Heritability:

     // create 4 LD-stratified GRMs based on stratifying SNPs into LD-score quartiles for calculating LD-adjusted heritability values
     CREATE_LDSTRATIFIED_GRMS(all_ldscore_results, FILTER_PLINKSET.out, params.autosome_num)

     // calculate LD-adjusted heritability (constrained)
     CALCULATE_HERITABILITY_REML_LDSTRATIFIED_CONSTRAINED(CREATE_LDSTRATIFIED_GRMS.out, GENERATE_PHENOFILE.out.flatten(), covar, qcovar)

    // calculate LD-adjusted heritability (no constraint)
     CALCULATE_HERITABILITY_REML_LDSTRATIFIED_NOCONSTRAINT(CREATE_LDSTRATIFIED_GRMS.out, GENERATE_PHENOFILE.out.flatten(), covar, qcovar)
    

    // calculate heritability (constrained)
    CALCULATE_HERITABILITY_REML_CONSTRAINED(CREATE_FULL_GRM.out, GENERATE_PHENOFILE.out.flatten(), covar, qcovar)

    // calculate heritability (no constraint)
    CALCULATE_HERITABILITY_REML_NOCONSTRAINT(CREATE_FULL_GRM.out, GENERATE_PHENOFILE.out.flatten(), covar, qcovar)
    
    
    // Genome-wide association studies (GWAS), clumping results, and plotting:

    // run GWAS MLMA LOCO (mixed linear model association leave-one-chromosome-out approach)
    GWAS_MLMA_LOCO(FILTER_PLINKSET.out, CREATE_FULL_GRM.out, gwas_chr_ch, params.autosome_num, GENERATE_PHENOFILE.out.flatten(), covar, qcovar)

    // create a tuple of per chr GWAS MLMA files grouped by phenotype
    gwas_results_tuple_ch = GWAS_MLMA_LOCO.out.groupTuple()
                                    
 
    // merge per chr GWAS results into a single GWAS output file for each phenotype analyzed
    MERGE_GWAS_RESULTS(gwas_results_tuple_ch)

   // create a tuple of GWAS .mlma.loco output file matched by phenotype 
   gwas_output_tuple_ch = MERGE_GWAS_RESULTS.out 

   // create a tuple of phenotype files to combine with tuple of GWAS output files matched by phenotype name to be submitted to clump process
   // pheno_tuple_ch will have corresponding phenotype file (key will be the full file name as splitting is occuring on the '.' which is also not present in terra variable name)
   pheno_tuple_ch = GENERATE_PHENOFILE.out.flatMap().map { file -> def key = file.name.toString().tokenize('.').get(0)
                                                           return tuple(key,file) }
   
   // combine tuple of phenotype files and GWAS .mlma.loco output file matched by phenotype (common key)
   gwas_pheno_pairing_ch = pheno_tuple_ch.combine(gwas_output_tuple_ch, by:0)

   // combine gwas_pheno_pairing_ch with filtered plinkset files to be submitted to CLUMP process for each phenotype entry
   gwas_pheno_pairing_plinkset_ch = gwas_pheno_pairing_ch.combine(FILTER_PLINKSET.out)


   // clump GWAS results 
   CLUMP(gwas_pheno_pairing_plinkset_ch, params.clump_pval, params.clump_kb, params.clump_rsquared, params.autosome_num, params.annotation)

   // save .clumped output file from clumping for submission to plotting process
   clumped_outfile_forplotting = CLUMP.out.flatten().filter { it -> it =~ /.clumped$/ }
    
   // clump channel with phenotype as key and corresponding clumped output files
   clump_tuple_ch = clumped_outfile_forplotting.map { file -> def key = file.name.toString().tokenize('.').get(0)
                                                      return tuple(key, file) }

   // pair clumped output file with GWAS output file based on phenotype name as key
   gwas_clump_pairing_ch = clump_tuple_ch.combine(gwas_output_tuple_ch, by:0)

   // plot manhattan and qqplots to display GWAS results for each phenotype of interest
   PLOT_MANHATTAN_QQPLOT(gwas_clump_pairing_ch) 

 
}

