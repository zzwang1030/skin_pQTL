#!/bin/bash
#SBATCH --nodes=1                   
#SBATCH --ntasks-per-node=1        
#SBATCH --mem-per-cpu=49g
#SBATCH --partition=g2_user        
#SBATCH --job-name=fusion_assoc_skin_protein_bslmm    
#SBATCH --output=skin_protein_fusion_assoc.bslmm.%A.%a.out        
#SBATCH --error=skin_protein_fusion_assoc.bslmm.%A.%a.err
#SBATCH --array=1-22

################# Software Envrironment ###########################
unset I_MPI_PMI_LIBRARY                        
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0   
export R_LIBS_USER="/sh2/home/sunyuanqiang/.conda/envs/python39/lib/R/library" 
module load intel/2022 intelmpi/2022          
module load R
# module load plink2

################ Run job ###################################
export OMP_NUM_THREADS=1
CHR=${SLURM_ARRAY_TASK_ID}
OUT="leprosy_GWAS_metadata_lwc_20240204.pos_id.fusion.skin_protein.${CHR}.bslmm.txt"

Rscript ~/software/fusion_twas-master/FUSION.assoc_test.syq.R --sumstats leprosy_GWAS_metadata_lwc_20240204.pos_id.ldsc.sumstats.gz --weights skin_protein.bslmm.wgt.pos --ref_ld_chr /sh2/home/sunyuanqiang/reference/G1000_hg38/g1000_eas_hg38.pos_id. --chr ${CHR} --out ${OUT}
