#!/bin/bash
#SBATCH -p shared                                               # Partition to submit to
#SBATCH -n 16                                           # Number of cores
#SBATCH -N 1                                            # Ensure that all cores are on one machine
#SBATCH -t 3-0:00                                       # Runtime in days-hours:minutes
#SBATCH --mem 100000                                    # Memory in MB
#SBATCH -J gstacks                                  # job name
#SBATCH -o gstacks%A.out                   # File to which standard out will be written
#SBATCH -e gstacks%A.err                   # File to which standard err will be written
#SBATCH --mail-type=ALL                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=schaturvedi@fas.harvard.edu     # Email to which notifications will be sent

# load modules
module purge
module load gcc/7.1.0-fasrc01 stacks/2.4-fasrc01

#specify file paths and directories
#specify paths

#path to the genome assembly, note you only need to specify the prefix
bwa_db=/n/holyscratch01/hopkins_lab/Chaturvedi/ddrad_sam_april2021/genome_assembly/curated.with.repeats

#path to write the output stack files
mkdir /n/holyscratch01/hopkins_lab/Chaturvedi/ddrad_sam_april2021/alignment_varcalling/stacks
stacks_out=/n/holyscratch01/hopkins_lab/Chaturvedi/ddrad_sam_april2021/alignment_varcalling/stacks

#path to bamfiles
bam_dir=/n/holyscratch01/hopkins_lab/Chaturvedi/ddrad_sam_april2021/alignment_varcalling/aligned_bamfiles


# Run gstacks to build loci from the aligned paired-end data. We have instructed
# gstacks to remove any PCR duplicates that it finds.

gstacks -I $bam_dir -M $bam_dir/popmap_geography.txt --rm-pcr-duplicates -O $stacks -t 8

# Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics and 
# smooth the statistics across the genome. Export several output files.

populations -P $stacks -M $bam_dir/popmap_geography.txt -r 0.65 --vcf --genepop --fstats --smooth --hwe -t 8
