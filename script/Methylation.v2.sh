#
/usr/bin/bash
#PBS -N Methylation.v2
#PBS -l walltime=120:00:00
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
SUFFIX=fastq.gz
rand=$RANDOM
FastqPath=/fs/project/PAS1475/Yuzhou_Chang/Methylation/Amer_Fastq/CleanFq

cd /fs/project/PAS1475/Yuzhou_Chang/Methylation/Amer_Fastq/New_analysis_Bowtie2
MYPATH=${pwd}
mkdir Sorted_Bam
# executive script
#mkdir CleanFq
ls ${FastqPath}/*fq.gz |cut -d/ -f9|cut -d_ -f1|uniq > mylist
# generate G-->A or C-->T genome fasta, this should be run for only one times. Once finish running it, please comment it. 
bismark_genome_preparation --verbose --bowtie2 --path_to_aligner /users/PAS1475/cyz931123/software/bowtie2-2.3.5.1-linux-x86_64 --parallel 10 /users/PAS1475/cyz931123/reference/reference
# do loop work for quality control and alignment.
for i in `cat mylist`
do
mkdir -p ${i}/QC_BeforeClean ${i}/QC_AfterClean ${i}/Results ${i}/Results/Extract ${i}/Results/Report Link2Bam
# quality control and check
fastqc -o ${i}/QC_BeforeClean -f fastq ${i}*R1.${SUFFIX} ${i}*R2.${SUFFIX}

${TRIM} --paired -q 20 --rrbs --illumina -o CleanFq --gzip --path_to_cutadapt /users/PAS1475/cyz931123/.local/bin/cutadapt ${i}*R1.${SUFFIX} ${i}*R2.${SUFFIX}

fastqc -o ${i}/QC_AfterClean -f fastq CleanFq/${i}*1.fq.gz CleanFq/${i}*2.fq.gz

mkdir -p ${i}/Results ${i}/Results/Extract ${i}/Results/Report Link2Bam
#alignment 
bismark -q --path_to_bowtie2 /users/PAS1475/cyz931123/software/bowtie2-2.3.5.1-linux-x86_64 --N 1 --L 15 -D 50 --score_min L,-0.6,-0.6 -p 10 --X 600 --gzip --bowtie2 --no-spliced-alignment -o ${i}/Results --unmapped --ambiguous --samtools_path /users/PAS1475/cyz931123/software/samtools-1.9 --genome_folder /users/PAS1475/cyz931123/reference/reference -1 ${FastqPath}/${i}*1.fq.gz -2 ${FastqPath}/${i}*2.fq.gz

# extract information from bismark output
bismark_methylation_extractor --paired-end --comprehensive -o ${i}/Results/Extract --samtools_path /users/PAS1475/cyz931123/software/samtools-1.9 --gzip --cytosine_report --genome_folder /users/PAS1475/cyz931123/reference/reference --parallel 10 --bedGraph --buffer_size 20G ${i}/Results/${i}*bam 
# report extract and quality check
bismark2report --dir ${i}/Results/Report --alignment_report ${i}/Results/${i}*_report.txt --splitting_report ${i}/Results/Extract/${i}*_splitting_report.txt --mbias_report ${i}/Results/Extract/${i}*M-bias.txt 
ln -s ${MYPATH}/${i}/Results/${i} ${MYPATH}/Link2Bam/
/users/PAS1475/cyz931123/software/samtools-1.9/samtools view -h ${i}/Results/${i}*.bam | /users/PAS1475/cyz931123/software/samtools-1.9/samtools view -Sh -| /users/PAS1475/cyz931123/software/samtools-1.9/samtools sort -T Sorted_BAM/${i}_${rand} -o Sorted_BAM/${i}.sort.bam -O bam -
done

# Globally report alignment 
bismark2summary ${MYPATH}/Amer_Fastq/Link2Bam/*bam
