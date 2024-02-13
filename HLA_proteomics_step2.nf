#!/usr/bin/env nextflow

params.protein_data = "/well/ckb/users/aey472/projects/HLA_proteomics/data"
params.genetic_data = "/well/ckb/users/aey472/projects/HLA_proteomics/data/genetic_data"
params.output = "/well/ckb/users/aey472/projects/HLA_proteomics/output"
params.bedfile = "/well/ckb/shared/filesystem/genetic_data/GWAS_data/b38_bpca/b38_bpca_100706" 

protein_ids = file("${params.protein_data}/GID_proteins.txt")

HLA_imputed = "/well/ckb-share/CKB_HLA_imputed_V2/combined/CKB_HLA_imputed_V2.newnames.GWAS_inds.vcf.gz"

Channel
    .fromPath("${params.protein_data}/olink_proteins_3000.txt")
    .splitCsv(header: false, sep: '\t')
    .map { row -> row[0] }
    .set { proteins }

process RegenieStep1 {
    tag "${protein}"
    publishDir path: params.output, mode: 'copy'

    input:
    val protein from proteins

    output:
    file "fit_bin_out_${protein}" into regenie_step1_output

    script:
    """
    mamba activate regenie
    regenie \\
        --step 1 \\
        --bed ${params.bedfile} \\
        --covarFile ${params.protein_data}/data_baseline_covariates_regenie.txt \\
        --phenoFile ${params.protein_data}/data_baseline_olink_allPanels_regenie.txt \\
        --phenoCol ${protein} \\
        --covarColList gwas_array_type,national_pc01,national_pc02,national_pc03,national_pc04,national_pc05,national_pc06,national_pc07,national_pc08,national_pc09,national_pc10,national_pc11 \\
        --catCovarList gwas_array_type \\
        --bsize 100 \\
        --qt \\
        --force-qt \\
        --out ${params.output}/fit_bin_out_${protein}
    """
}

process RegenieStep2 {
    tag "${protein}"
    publishDir path: params.output, mode: 'copy'

    input:
    val protein from regenie_step1_output

    script:
    """
    mamba activate regenie
    regenie \\
        --step 2 \\
        --pgen ${params.genetic_data}/CKB_HLA_imputed_V2.newnames.GWAS_inds.protein_inds \\
        --covarFile ${params.protein_data}/data_baseline_covariates_regenie.txt \\
        --phenoFile ${params.protein_data}/data_baseline_olink_allPanels_regenie.txt \\
        --bsize 200 \\
        --covarColList gwas_array_type,national_pc01,national_pc02,national_pc03,national_pc04,national_pc05,national_pc06,national_pc07,national_pc08,national_pc09,national_pc10,national_pc11 \\
        --catCovarList gwas_array_type \\
        --phenoCol ol_pdia4 \\
        --qt \\
        --firth \\
        --pThresh 0.01 \\
        --pred ${params.output}/fit_bin_out_${protein}.list \\
        --out ${params.output}/test_bin_out_${protein}_firth
    """
}
