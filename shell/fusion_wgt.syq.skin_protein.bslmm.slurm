#!/bin/bash
#SBATCH --nodes=1                   # 节点数量
#SBATCH --ntasks-per-node=1        # 每个节点核心数量，每个节点最多 56 个核心
#SBATCH --mem-per-cpu=20g
#SBATCH --partition=g2_user        # 队列分区且必须指定正确分区
#SBATCH --job-name=fusion_wgt_skin_protein_bslmm    # 作业名称
#SBATCH --output=skin_protein_fusion_wgt_bslmm.%A.%a.out        # 正常日志输出 (%j 参数值为 jobId)
#SBATCH --error=skin_protein_fusion_wgt_bslmm.%A.%a.err         # 错误日志输出 (%j 参数值为 jobId)
#SBATCH --array=1-52

################# Software Envrironment ###########################
unset I_MPI_PMI_LIBRARY                        # 取消默认mpi库，使用intel自带
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0   # intel 多节点作业所需修改参数
module load intel/2022 intelmpi/2022           # intel 环境加载
# module load plink2

################ Run job ###################################
export OMP_NUM_THREADS=1
bash fusion_wgt.syq.new.bslmm.sh skin_protein skin_protein_all_qqnorm.txt.gz skin_protein_all_qqnorm.txt.gz.HEADER skin_protein_all_qqnorm.txt.gz.12pheno_5geno.PCs.plink.covar /sh2/home/sunyuanqiang/projects/pQTL/formal/hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.plink
