#!/bin/bash
#SBATCH -A ckb.prj
#SBATCH -J HLA_proteomics
#SBATCH -o HLA_proteomics_%A_%a.out
#SBATCH -e HLA_proteomcis_%A_%a.err
#SBATCH -p short
#SBATCH -c 2
#SBATCH -a 1-2930


protein_data=/well/ckb/users/aey472/projects/HLA_proteomics/data
protein_ids=${protein_data}/GID_proteins.txt

HLA_imputed=/well/ckb-share/CKB_HLA_imputed_V2/combined/CKB_HLA_imputed_V2.newnames.GWAS_inds.vcf.gz

data=/well/ckb/users/aey472/projects/HLA_proteomics/data
genetic_data=/well/ckb/users/aey472/projects/HLA_proteomics/data/genetic_data

output=/well/ckb/users/aey472/projects/HLA_proteomics/output

protein=$(sed -n ${SLURM_ARRAY_TASK_ID}'{p;q}' ${data}/olink_proteins_3000.txt)

mamba activate regenie

regenie \
        --step 1 \
        --bed ${bedfile} \
        --covarFile ${data}/data_baseline_covariates_regenie.txt \
        --phenoFile ${data}/data_baseline_olink_allPanels_regenie.txt \
				--phenoCol ${protein} \
        --covarColList gwas_array_type,national_pc01,national_pc02,national_pc03,national_pc04,national_pc05,national_pc06,national_pc07,national_pc08,national_pc09,national_pc10,national_pc11 \
        --catCovarList gwas_array_type \
        --bsize 100 \
        --qt \
        --force-qt \
        --out ${output}/fit_bin_out_${protein}

regenie \
        --step 2 \
        --pgen ${genetic_data}/CKB_HLA_imputed_V2.newnames.GWAS_inds.protein_inds \
        --covarFile ${data}/data_baseline_covariates_regenie.txt \
        --phenoFile ${data}/data_baseline_olink_allPanels_regenie.txt \
        --bsize 200 \
        --covarColList gwas_array_type,national_pc01,national_pc02,national_pc03,national_pc04,national_pc05,national_pc06,national_pc07,national_pc08,national_pc09,national_pc10,national_pc11 \
        --catCovarList gwas_array_type \
        --phenoCol ol_pdia4 \
        --qt \
        --firth \
        --pThresh 0.01 \
        --pred ${output}/fit_bin_out_${protein}.list \
        --out ${output}/test_bin_out_${protein}_firth
