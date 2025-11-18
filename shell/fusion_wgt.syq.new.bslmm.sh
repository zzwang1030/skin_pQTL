#!/bin/bash
set -e

# 检查参数数量
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <PRE> <NormExpression> <ExpressionHeader> <COVAR> <GENOTYPE>"
    exit 1
fi

# 从命令行参数读取变量
PRE=$1
NormExpression=$2
ExpressionHeader=$3
COVAR=$4
GENOTYPE=$5
NR=${SLURM_ARRAY_TASK_ID}

# 工具路径
GCTA="$HOME/software/fusion_twas-master/gcta_nr_robust"
PLINK="plink --allow-no-sex"
GEMMA="$HOME/software/gemma-0.98.5-linux-static/gemma-0.98.5-linux-static-AMD64"
FUSION="$HOME/software/fusion_twas-master"

# Input files needed:
# zcat $NormExpression | head -n1 | tr '\t' '\n' > $NormExpression.HEADER # 这个会输出到对应的绝对路径下
# datamash transpose < $COVARaw > mytmp
# awk '{print $1"\t"$0}' mytmp > mytmp2
# awk 'NR==1{$1="FID"; $2="IID"}1' OFS="\t" mytmp2 > $COVARaw.plink.covar # 这个会输出到对应的绝对路径下
# COVAR=$COVARaw.plink.covar

mkdir -p tmp_bslmm
mkdir -p WEIGHTS_bslmm
mkdir -p HSQ_bslmm

zcat $NormExpression | tail -n+2 | awk -v i=$NR 'NR > (i-1)*100 && NR <= i*100' |  while read PARAM; do

    # PARAM contains the expression values for this gene
    GNAME=$(echo $PARAM | awk '{ print $4 }')
    CHR=$(echo $PARAM | awk '{ print $1 }')
    P0=$(echo $PARAM | awk '{ p=$2 - 1000e3; if(p<0) p=0; print p; }')
    P1=$(echo $PARAM | awk '{ print $3 + 1000e3 }')
    echo $GNAME $CHR $P0 $P1

    # OUT="tmp/$PRE.$NR/$PRE.$GNAME"
    OUT="tmp_bslmm/$PRE.$GNAME"
    # Convert the expression matrix to a PLINK format phenotype
    echo $PARAM | tr ' ' '\n' | paste $ExpressionHeader $ExpressionHeader - | tail -n+5 > $OUT.pheno || continue
    # Extract the locus around this gene
    $PLINK --silent --bfile $GENOTYPE --chr $CHR --from-bp $P0 --to-bp $P1 --make-bed --out $OUT --pheno $OUT.pheno --keep $OUT.pheno --maf 0.0001 || continue
    if [ ! -f $OUT.bed ]; then
        continue
    fi

    # Run FUSION
    FINAL_OUT="WEIGHTS_bslmm/$PRE.$GNAME"
    Rscript $FUSION/FUSION.compute_weights.R --bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT --verbose 0 --save_hsq --PATH_gcta $GCTA --PATH_gemma $GEMMA --models blup,lasso,top1,enet,bslmm --covar $COVAR  --hsq_p 1.0 --crossval 5 || continue # used bslmm model, but slowly

    cat $FINAL_OUT.hsq | awk -v w="$PRE $GNAME" '{ print w,$0 }' > HSQ_bslmm/$PRE.$GNAME.hsq
    echo "$(pwd)/$FINAL_OUT.wgt.RDat $GNAME $CHR $P0 $P1" | tr ' ' '\t'  >> $PRE.bslmm.wgt.pos # this output pos file is needed in the next associated step. by syq
    
    rm -f $FINAL_OUT.hsq
    rm $OUT.*
done
