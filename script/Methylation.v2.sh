#!/usr/bin/bash
#PBS -N Methylation.v2
#PBS -l walltime=80:00:00
#PBS -l nodes=1:ppn=28
#PBS -l mem=64GB
#PBS -m abe
#PBS -j oe
#PBS -S /usr/bin/bash
set -x

#load mofule
module load fastqc

# set global environment 
export PATH=$PATH:/users/PAS1475/cyz931123/software/bismark_v0.22.1
TRIM=/users/PAS1475/cyz931123/software/TrimGalore-0.6.1/trim_galore
MYPATH=$(pwd)
SUFFIX=fastq.gz
rand=$RANDOM

cd Amer_Fastq
# mkdir Sorted_BAM
# executive script
# mkdir CleanFq
# ls *.fastq.gz |cut -d_ -f1|uniq > mylist
# generate G-->A or C-->T genome fasta, this should be run for only one times. Once finish running it, please comment it. 
#bismark_genome_preparation --verbose --hisat2 --path_to_aligner /BigData/software/hisat2-2.1.0 --parallel 4 /BigData/reference/human
# do loop work for quality control and alignment.
for i in `cat mylist2`
do
# mkdir -p ${i}/QC_BeforeClean ${i}/QC_AfterClean ${i}/Results ${i}/Results/Extract ${i}/Results/Report Link2Bam
# quality control and check
# fastqc -o ${i}/QC_BeforeClean -f fastq ${i}*R1.${SUFFIX} ${i}*R2.${SUFFIX}

# ${TRIM} --paired -q 20 --rrbs --illumina -o CleanFq --gzip --path_to_cutadapt /users/PAS1475/cyz931123/.local/bin/cutadapt ${i}*R1.${SUFFIX} ${i}*R2.${SUFFIX}

# fastqc -o ${i}/QC_AfterClean -f fastq CleanFq/${i}*1.fq.gz CleanFq/${i}*2.fq.gz
# fix the bug
rm -rf ${i}/Results
mkdir -p ${i}/Results ${i}/Results/Extract ${i}/Results/Report Link2Bam
#alignment 
bismark -q --path_to_hisat2 /users/PAS1475/cyz931123/software/hisat2-2.1.0 --gzip --hisat2 --no-spliced-alignment -o ${i}/Results --unmapped --ambiguous --samtools_path /users/PAS1475/cyz931123/software/samtools-1.9 /users/PAS1475/cyz931123/reference/reference -1 CleanFq/${i}*1.fq.gz -2 CleanFq/${i}*2.fq.gz

# extract information from bismark output
bismark_methylation_extractor --paired-end --comprehensive -o ${i}/Results/Extract --samtools_path /users/PAS1475/cyz931123/software/samtools-1.9 --gzip --cytosine_report --genome_folder /users/PAS1475/cyz931123/reference/reference --parallel 10 --bedGraph --buffer_size 20G ${i}/Results/${i}*bam 
# report extract and quality check
bismark2report --dir ${i}/Results/Report --alignment_report ${i}/Results/${i}*_report.txt --splitting_report ${i}/Results/Extract/${i}*_splitting_report.txt --mbias_report ${i}/Results/Extract/${i}*M-bias.txt 
ln -s ${MYPATH}/Amer_Fastq/${i}/Results/${i}* ${MYPATH}/Amer_Fastq/Link2Bam/
/users/PAS1475/cyz931123/software/samtools-1.9/samtools view -h ${i}/Results/${i}*.bam | /users/PAS1475/cyz931123/software/samtools-1.9/samtools view -Sh -| /users/PAS1475/cyz931123/software/samtools-1.9/samtools sort -T Sorted_BAM/${i}_${rand} -o Sorted_BAM/${i}.sort.bam -O bam -
done

# Globally report alignment 
bismark2summary ${MYPATH}/Amer_Fastq/Link2Bam/*bam




