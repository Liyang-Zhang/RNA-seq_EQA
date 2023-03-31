#!/bin/bash

#SBATCH --mail-type=end
#SBATCH --mail-user=zly18810602991@163.com
#SBATCH --job-name=test
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -J rna
#SBATCH --ntasks-per-node=8
#SBATCH --output=%j.o
#SBATCH --error=%j.e

# Initialize a variable with an intuitive name to store the name of the input fastq file
R1=$1
R2=$2
sample=$3
project_dir=$4

# Data from Pfizer
#baseSuper=`basename ${sample} .fastq`
#baseSuper=`basename ${sample} .fastq.gz`
baseSuper=$(echo ${sample} | cut -d'-' -f 3)
base=${baseSuper:3}
#base=$(echo ${baseSuper} | cut -d'_' -f 1)  # get the filename with basename function
echo ${baseSuper}

# specify the number of cores to use
cores=8

# reference files
#transcriptSeq=/dssg/home/acct-medkwf/medkwf4/ngs/hbctraining/rnaseq/reference/Homo_sapiens.GRCh38.cdna.all.fa
#chr1_dna=/dssg/home/acct-medkwf/medkwf4/ngs/hbctraining/rnaseq/reference/Homo_sapiens.GRCh38.dna.chromosome.1.fa
gtf=/dssg/home/acct-medkwf/medkwf4/reference/rna_EQA/Homo_sapiens.GRCh38.106.spikein.gtf
dna=/dssg/home/acct-medkwf/medkwf4/reference/rna_EQA/Homo_sapiens.GRCh38.dna.primary_assembly.spikein.fa
transcriptSeq=/dssg/home/acct-medkwf/medkwf4/reference/rna_EQA/Homo_sapiens.GRCh38.cdna.all.spikein.fa

# set up output filenames and locations
fastqc_out=${project_dir}/results/fastqc
trim_galore_out=${project_dir}/results/trim_galore
#bams_out=${project_dir}/results/STAR/bams
bams_out=${project_dir}/bams
#align_out=${project_dir}/results/STAR/${sample}_
align_out=${project_dir}/results/STAR/${base}
sample_bam=${project_dir}/results/STAR/${base}Aligned.sortedByCoord.out.bam
#qualimap_out=${project_dir}/results/qualimap
qualimap_out=${project_dir}/qualimap
counts=${project_dir}/results/counts/${base}_featurecounts.txt
counts_out=${project_dir}/results/counts
stringtie_out=${project_dir}/results/stringtie
multiqc=${project_dir}/multiqc
salmon_out=${project_dir}/results/salmon

# make all of the output directories
mkdir -p ${project_dir}/results/fastqc
mkdir -p ${project_dir}/results/trim_galore
mkdir -p ${project_dir}/results/STAR/bams
mkdir -p ${project_dir}/bams
#mkdir -p ${project_dir}/results/qualimap
mkdir -p ${project_dir}/qualimap
mkdir -p ${project_dir}/results/counts
mkdir -p ${project_dir}/results/stringtie
mkdir -p ${multiqc}
mkdir -p ${salmon_out}

# make and specify directories 
STAR_index=/dssg/home/acct-medkwf/medkwf4/reference/rna_EQA/STAR_hg38_index
salmon_index=/dssg/home/acct-medkwf/medkwf4/reference/rna_EQA/salmon_hg38_index

# set up the environment
module load miniconda3
source activate ngs-rna

# Run FastQC for raw
#fastqc -t ${cores} \
#${R1} \
#-o ${fastqc_out}
#
#fastqc -t ${cores} \
#${R2} \
#-o ${fastqc_out}

# trim adapters
#trim_galore --quality 20 --phred33 --fastqc --length 20 --gzip \
#--cores 4 \
#-o ${trim_galore_out} \
#--paired ${R1} ${R2} \
#2> ${trim_galore_out}/${sample}.trimgalore.log

# use salmon for mapping and quantification
# INDEX ONLY ONCE
#salmon index \
#-t ${transcriptSeq} \
#-i ${salmon_index} \
#-k 31

salmon quant -i ${salmon_index} \
-l A \
-1 ${trim_galore_out}/${base}.trimmed.R1.fastq.gz \
-2 ${trim_galore_out}/${base}.trimmed.R2.fastq.gz \
-o ${salmon_out}/${base}.salmon \
--useVBOpt \
--seqBias \
--validateMappings \
-p ${cores}

# use STAR to make index and mapping
# INDEX ONLY ONCE
#STAR --runThreadN 6 \
#--runMode genomeGenerate \
#--genomeDir ${STAR_index} \
#--genomeFastaFiles ${dna} \
#--sjdbGTFfile ${gtf} \
#--sjdbOverhang 99 \
#--genomeSAindexNbases 12 

#mv ${trim_galore_out}/${sample}_combined_R1_val_1.fq.gz ${trim_galore_out}/${base}.trimmed.R1.fastq.gz
#mv ${trim_galore_out}/${sample}_combined_R2_val_2.fq.gz ${trim_galore_out}/${base}.trimmed.R2.fastq.gz
#mv ${trim_galore_out}/${sample}_combined_R1_val_1_fastqc.html ${trim_galore_out}/${base}.trimmed.R1.fastqc.html
#mv ${trim_galore_out}/${sample}_combined_R2_val_2_fastqc.html ${trim_galore_out}/${base}.trimmed.R2.fastqc.html
#mv ${trim_galore_out}/${sample}_combined_R1_val_1_fastqc.zip ${trim_galore_out}/${base}.trimmed.R1.fastqc.zip
#mv ${trim_galore_out}/${sample}_combined_R2_val_2_fastqc.zip ${trim_galore_out}/${base}.trimmed.R2.fastqc.zip

#STAR --genomeDir ${STAR_index} \
#--runThreadN ${cores} \
#--readFilesIn ${trim_galore_out}/${base}.trimmed.R1.fastq.gz ${trim_galore_out}/${base}.trimmed.R2.fastq.gz \
#--readFilesCommand zcat \
#--outFileNamePrefix ${align_out} \
#--outSAMtype BAM SortedByCoordinate \
#--outSAMunmapped Within \
#--outSAMattributes Standard 

# Move ${sample}*.out.bam to ${bams_out} or create a softlink for bam files
#mv ${sample_bam} ${bams_out}
#ln -s ${sample_bam} ${bams_out}

# rename and create BAM index as required
#mv ${bams_out}/${base}Aligned.sortedByCoord.out.bam ${bams_out}/${base}.bam
#samtools index ${bams_out}/${base}.bam

# Obtain quality metrics using Qualimap
#EQA sequencing is non-strand-specific protocol, so -p should not be set
#-p strand-specific-reverse
#unset DISPLAY
#
#qualimap rnaseq \
#-outdir ${qualimap_out}/${base} \
#-a proportional \
#-bam ${bams_out}/${base}.bam \
#-gtf ${gtf} \
#--java-mem-size=24G \
#2> ${qualimap_out}/${base}_qualimap.screenoutput


#multiqc -n multiqc_report_rnaseq_EQA \
#-o ${multiqc} \
#${trim_galore_out} \
#${project_dir}/results/STAR/*Log.final.out \
#${qualimap_out} \
#/dssg/home/acct-medkwf/medkwf4/ngs/hbctraining/rnaseq/results/salmon/*

# Count mapped reads in two ways as #1 and #2 below
#1 Create read count for the single file of $sample
#featureCounts -T ${cores} -s 1 \
#-a ${gtf} \
#-o ${counts} \
#${bams_out}/${base}Aligned.sortedByCoord.out.bam \
#2> ${counts_out}/${base}_featurecounts.screenoutput

#featureCounts -T ${cores} -p -s 0 \
#-a ${gtf} \
#-o ${counts} \
#${bams_out}/${base}.bam \
#2> ${counts_out}/${base}_featurecounts.screenoutput


#2 Create read counts for all bam files in the directory of ${bams_out}. Comment this command until you get bam files for all 16 data.

#featureCounts -T ${cores} -p -s 0 \
#-a ${gtf} \
#-o ${counts_out}/all_featurecounts_protein_coding.txt \
#${bams_out}/*.bam \
#2> ${counts_out}/all_featurecounts.screenoutput

# Stringtie for FPKM, TPM
#--rf   Assumes a stranded library fr-firststrand
#--fr   Assumes a stranded library fr-secondstrand
#stringtie -v -p ${cores} -e --fr \

#stringtie -v -p ${cores} -e \
#-o ${stringtie_out}/${base}.gtf -G ${gtf} \
#-A ${stringtie_out}/${base}.gene.abund.txt ${bams_out}/${base}.bam

#stringtie -v -p ${cores} -e \
#-o ${stringtie_out}/${sample}.gtf -G ${gtf} \
#-A ${stringtie_out}/${sample}.gene.abund.txt ${bams_out}/${sample}*Aligned.sortedByCoord.out.bam


#Transcript merge mode
#Not we want
#mkdir -p ${stringtie_out}/merge
#stringtie --merge \
#-o ${stringtie_out}/merge/EQA_transcript.gtf \
#-p ${cores} \
#-G ${gtf} \

#mkdir -p ${stringtie_out}/merge
#stringtie_merge_transcript=${stringtie_out}/merge/EQA_transcript.txt
#stringtie_merge_gene=${stringtie_out}/merge/EQA_gene.txt
#count_transcript=${stringtie_out}/merge/count_transcript.txt
#count_gene=${stringtie_out}/merge/count_gene.txt
#rm -fr ${stringtie_merge_transcript} ${stringtie_merge_gene}
#touch ${stringtie_merge_transcript}
#touch ${stringtie_merge_gene}
#
#for gtf_file in `ls -l ${stringtie_out} | grep .gtf | awk '{print  $9}' | uniq`
#do
#  echo "GTF file is ${stringtie_out}/${gtf_file}"
#  sample_name=$(echo ${gtf_file} | cut -d'.' -f 1)
#  grep -P "\ttranscript\t" ${stringtie_out}/${gtf_file} | cut -f9 | awk -v sample=$sample_name '{gsub("\"","",$0);gsub(";","",$0);print sample,$4,$(NF-2),$NF}' >> ${stringtie_merge_transcript}
#done
#
#for gene_file in `ls -l ${stringtie_out} | grep .gene.abund.txt | awk '{print  $9}' | uniq`
#do
#  echo "gene abundance file is ${stringtie_out}/${gene_file}"
#  sample_name=$(echo ${gene_file} | cut -d'.' -f 1)
#  grep -vE "^Gene" ${stringtie_out}/${gene_file} | awk -v sample=$sample_name '{print sample,$1,$(NF-1),$NF}' >> ${stringtie_merge_gene}
#done

# Use the official python3 script to get the gene and transcript level counts table
# the stringtie_lst.txt lists the name and path of the stringtie output GTF files
#python3 /dssg/home/acct-medkwf/medkwf4/script/rna/prepDE.py3 -i /dssg/home/acct-medkwf/medkwf4/script/rna/stringtie_lst.txt -t ${count_transcript} -g ${count_gene} 


