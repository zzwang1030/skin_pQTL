#!/bin/bash
#SBATCH --nodes=1                   # 节点数量
#SBATCH --ntasks-per-node=1        # 每个节点核心数量，每个节点最多 56 个核心
#SBATCH --mem-per-cpu=49g
#SBATCH --partition=g2_user        # 队列分区且必须指定正确分区
#SBATCH --job-name=fusion_assoc_skin_protein_bslmm    # 作业名称
#SBATCH --output=skin_protein_fusion_assoc.bslmm.%A.%a.out        # 正常日志输出 (%j 参数值为 jobId)
#SBATCH --error=skin_protein_fusion_assoc.bslmm.%A.%a.err         # 错误日志输出 (%j 参数值为 jobId)
#SBATCH --array=1-22

################# Software Envrironment ###########################
unset I_MPI_PMI_LIBRARY                        # 取消默认mpi库，使用intel自带
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0   # intel 多节点作业所需修改参数
export R_LIBS_USER="/sh2/home/sunyuanqiang/.conda/envs/python39/lib/R/library" # syq add: 防止单独跑不报错，sbatch 运行报错：不存在叫'plink2R'这个名字的程序包
module load intel/2022 intelmpi/2022           # intel 环境加载
module load R
# module load plink2

################ Run job ###################################
export OMP_NUM_THREADS=1
CHR=${SLURM_ARRAY_TASK_ID}
OUT="leprosy_GWAS_metadata_lwc_20240204.pos_id.fusion.skin_protein.${CHR}.bslmm.txt"

Rscript ~/software/fusion_twas-master/FUSION.assoc_test.R --sumstats leprosy_GWAS_metadata_lwc_20240204.pos_id.ldsc.sumstats.gz --weights skin_protein.bslmm.wgt.pos --ref_ld_chr /sh2/home/sunyuanqiang/reference/G1000_hg38/g1000_eas_hg38.pos_id. --chr ${CHR} --out ${OUT}
