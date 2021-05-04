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
cd /n/holyscratch01/hopkins_lab/Chaturvedi/ddrad_sam_april2021/alignment_varcalling/aligned_bamfiles/

for i in *.bam; do samtools idxstats $i | awk '{s+=$3} END {print s}'; done

#the slurm *out file is where all the numbers will be written
