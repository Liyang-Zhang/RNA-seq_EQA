#!/bin/bash
# raw data has been instored at the following  ${project_dir} already
# USAGE : bash $0 <project directory>
# hpc training data project_dir=/dssg/home/acct-medkwf/medkwf4/data/ncbi/public/sra/acclist/GSE50499/raw
# EQA: /dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/CG0701-R22024472-202208251343-1
# EQA cn500: /dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/cn500
project_dir=$1

#for sample in `ls -l ${project_dir} | awk '{print  $9}' | grep -E ".fastq.gz$" | uniq`
for sample in `ls -l ${project_dir} | awk '{print  $9}' | grep -E "R1.fastq.gz$" | uniq`
do
  base=$(echo ${sample} | cut -d'_' -f 1)
  R1=${project_dir}/${base}_R1.fastq.gz
  R2=${project_dir}/${base}_R2.fastq.gz
  #echo -e "Sample is ${sample/_1.fastq.gz/}"
  echo -e "Sample is ${base}"
 if
  [[ $sample == *"Lib202201"* ]]; then
    echo "running analysis on ${base}"
    sbatch /dssg/home/acct-medkwf/medkwf4/script/rna/pe_rnaseq.cn500.slurm ${R1} ${R2} ${base} ${project_dir}
  fi
done

