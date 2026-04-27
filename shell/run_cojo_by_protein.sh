#!/bin/bash

usage() {
  echo "Usage: $0 -i <input_file.gz> -b <bfile_prefix> [-d <diff_freq>] [-c <cojo_collinear>] [-p <cojo_p>] "
  echo "  -i input summary stats gz file (required)"
  echo "  -b plink bfile prefix (required)"
  echo "  -d --diff-freq parameter for COJO (default 0.4)"
  echo "  -c --cojo-collinear parameter for COJO (default 0.9)"
  echo "  -p --cojo-p parameter for COJO (default 7.90e-6)"
  exit 1
}


diff_freq=0.4
cojo_collinear=0.9
cojo_p=7.90e-6


while getopts ":i:b:d:c:p:" opt; do
  case $opt in
    i) INPUT="$OPTARG" ;;
    b) BFILE="$OPTARG" ;;
    d) diff_freq="$OPTARG" ;;
    c) cojo_collinear="$OPTARG" ;;
    p) cojo_p="$OPTARG" ;;
    *) usage ;;
  esac
done


if [[ -z "$INPUT" ]] || [[ -z "$BFILE" ]]; then
  echo "Error: input file and bfile prefix are required."
  usage
fi

mkdir -p tmp_cojo_input tmp_cojo_output

proteins=$(zcat "$INPUT" | tail -n +2 | cut -f1 | sort | uniq)

echo "Start processing proteins: $proteins"

for prot in $proteins; do
  echo "Processing protein: $prot"

  zcat "$INPUT" | awk -v prot="$prot" 'NR==1 || $1==prot' > tmp_cojo_input/${prot}.txt

  awk 'BEGIN{OFS="\t"; print "SNP","A1","A2","freq","b","se","p","N"} NR>1{
  split($2,a,":");
  print $2,a[4],a[3],$6,$8,$9,$7,189
}' tmp_cojo_input/${prot}.txt > tmp_cojo_input/${prot}_forCOJO.txt

  ~/software/gcta/gcta64 --bfile "$BFILE" \
    --cojo-file tmp_cojo_input/${prot}_forCOJO.txt \
    --cojo-slct \
    --diff-freq "$diff_freq" \
    --cojo-collinear "$cojo_collinear" \
    --cojo-p "$cojo_p" \
    --out tmp_cojo_output/${prot}_forCOJO_res
done

echo "All done."

