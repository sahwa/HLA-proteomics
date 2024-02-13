#!/usr/bin/env nextflow

params.protein_data = '/well/ckb/users/aey472/projects/HLA_proteomics/data'
params.protein_ids = "${params.protein_data}/GID_proteins.txt"
params.HLA_imputed = '/well/ckb-share/CKB_HLA_imputed_V2/combined/CKB_HLA_imputed_V2.newnames.GWAS_inds.vcf.gz'
params.data = '/well/ckb/users/aey472/projects/HLA_proteomics/data'
params.genetic_data = '/well/ckb/users/aey472/projects/HLA_proteomics/data/genetic_data'
params.output = '/well/ckb/users/aey472/projects/HLA_proteomics/output'

process HLAProteomics {
    tag "HLA_proteomics"
    
    publishDir "${params.output}", mode: 'copy'

    input:
    val protein_ids from Channel.fromPath(params.protein_ids)
    val HLA_imputed from Channel.value(params.HLA_imputed)
    val genetic_data from Channel.value(params.genetic_data)

    output:
    file("${genetic_data}/CKB_HLA_imputed_V2.newnames.GWAS_inds.protein_inds.*") into output_files

    script:
    """
    module purge all && module load BCFtools/1.17-GCC-12.2.0
    bcftools view -S \${protein_ids} \${HLA_imputed} -Oz > \${genetic_data}/CKB_HLA_imputed_V2.newnames.GWAS_inds.protein_inds.vcf.gz
    
    module purge all && module load PLINK/2.00a3.1-GCC-11.2.0
    plink2 --vcf \${genetic_data}/CKB_HLA_imputed_V2.newnames.GWAS_inds.protein_inds.vcf.gz --make-pgen --out \${genetic_data}/CKB_HLA_imputed_V2.newnames.GWAS_inds.protein_inds
    """
}

output_files.view { file -> 
    println("Generated file: ${file}")
}
