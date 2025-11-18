#!/bin/bash
#SBATCH --nodes=1                   
#SBATCH --ntasks-per-node=1        
#SBATCH --mem-per-cpu=20g
#SBATCH --partition=g2_user        
#SBATCH --job-name=fusion_wgt_skin_protein_bslmm    
#SBATCH --output=skin_protein_fusion_wgt_bslmm.%A.%a.out        
#SBATCH --error=skin_protein_fusion_wgt_bslmm.%A.%a.err         
#SBATCH --array=1-52

################# Software Envrironment ###########################
unset I_MPI_PMI_LIBRARY                       
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
module load intel/2022 intelmpi/2022
# module load plink2

################ Run job ###################################
export OMP_NUM_THREADS=1
bash fusion_wgt.syq.new.bslmm.sh skin_protein skin_protein_all_qqnorm.txt.gz skin_protein_all_qqnorm.txt.gz.HEADER skin_protein_all_qqnorm.txt.gz.12pheno_5geno.PCs.plink.covar /sh2/home/sunyuanqiang/projects/pQTL/formal/hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.plink
