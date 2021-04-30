#!/bin/bash
#SBATCH -p shared                                               # Partition to submit to
#SBATCH -n 16                                           # Number of cores
#SBATCH -N 1                                            # Ensure that all cores are on one machine
#SBATCH -t 3-0:00                                       # Runtime in days-hours:minutes
#SBATCH --mem 100000                                    # Memory in MB
#SBATCH -J bwaindex                                  # job name
#SBATCH -o bwaindex_%A.out                   # File to which standard out will be written
#SBATCH -e bwaindex_%A.err                   # File to which standard err will be written
#SBATCH --mail-type=ALL                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=schaturvedi@fas.harvard.edu     # Email to which notifications will be sent

# load a java module
module load bwa

#go to the directory with the genome assembly
cd /n/holyscratch01/hopkins_lab/Chaturvedi/ddrad_sam_april2021/genome_assembly/

#run bwa
bwa index -a bwtsw curated.with.repeats.fasta curated.with.repeats
