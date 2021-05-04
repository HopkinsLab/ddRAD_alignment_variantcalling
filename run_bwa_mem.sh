#!/bin/bash
#SBATCH -p shared                                               # Partition to submit to
#SBATCH -n 16                                           # Number of cores
#SBATCH -N 1                                            # Ensure that all cores are on one machine
#SBATCH -t 5-0:00                                       # Runtime in days-hours:minutes
#SBATCH --mem 150000                                    # Memory in MB
#SBATCH -J bwamem                                  # job name
#SBATCH -o bwamem_%A.out                   # File to which standard out will be written
#SBATCH -e bwamem_%A.err                   # File to which standard err will be written
#SBATCH --mail-type=ALL                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=schaturvedi@fas.harvard.edu     # Email to which notifications will be sent

# load modules
module load bwa
module load samtools
module load bcftools

#specify paths
#go to directory for fastq files
cd /n/holyscratch01/hopkins_lab/Chaturvedi/ddrad_sam_april2021/trimmomatic/trim_paired/
#path to the genome assembly, note you only need to specify the prefix
bwa_db=/n/holyscratch01/hopkins_lab/Chaturvedi/ddrad_sam_april2021/genome_assembly/curated.with.repeats.fasta
#path to write the output bam files
mkdir /n/holyscratch01/hopkins_lab/Chaturvedi/ddrad_sam_april2021/alignment_varcalling/aligned_bamfiles
out_path=/n/holyscratch01/hopkins_lab/Chaturvedi/ddrad_sam_april2021/alignment_varcalling/aligned_bamfiles

# Align paired-end data with BWA, convert to BAM and SORT.
for i in $(ls *.fastq | sed -r 's/_R[12]_paired.fastq//' | uniq)
do
        bwa mem -t 12 -k 15 -r 1.3 -T 20 $bwa_db  ${i}_R1_paired.fastq ${i}_R2_paired.fastq | samtools view -b | samtools sort --threads 10 > ${i}.bam
        samtools index ${i}.bam
	mv ${i}*bam* $out_path
done
