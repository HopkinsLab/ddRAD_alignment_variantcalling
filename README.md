## ddRAD_alignment_variantcalling_STACKS

This is the pipeline for alignment and variant calling using STACKS for Phlox ddRAD sequencing data.

**STEP 1: Alignmnet using bwa mem**
We will use BWA MEM algorithm to align the ddRAD data to the reference genome assembly. The first step of using BWA is to make an index of the reference genome in fasta format. The basic script for running this is (run_bwa_index.sh):

```bash
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
```
1. In the command above, I am specifying the 'bwstw' option for the index which is for bowtie. 
2. I maintain the prefix from the assembly name to make things easier for downstream processing.

The second step is to perform the actual alignment. The bwa mem algorithm is one of the three algorithms provided by BWA. It performs local alignment and produces alignments for different part of the query sequence. Here is the bash script for this alignment (run_bwa_mem.sh):


```bash
#!/bin/bash
#SBATCH -p shared                                               # Partition to submit to
#SBATCH -n 16                                           # Number of cores
#SBATCH -N 1                                            # Ensure that all cores are on one machine
#SBATCH -t 3-0:00                                       # Runtime in days-hours:minutes
#SBATCH --mem 100000                                    # Memory in MB
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
	mv ${i}.bam $out_path
done
```

**STEP 2: Run reference based alignment using STACKS**

First, we will create the population map file for using STACKS for alignment. This will be done on command line. We are creating a map based on geography as we are using the population IDs for downstream analyses. Read the STACKS manual to modify this for other factors such as phenotypes.

```bash
cd /n/holyscratch01/hopkins_lab/Chaturvedi/ddrad_sam_april2021/alignment_varcalling/aligned_bamfiles

# 1. Create list of samples ids for the aligned bams, we just need the prefix
ls *.fastq | sed -r 's/_R[12]_paired.fastq//' | uniq > sampleids.txt

# 2. Create list of populations
ls *.fastq | sed -r 's/_R[12]_paired.fastq//' | uniq | cut -d"-" -f 1 > populations.txt

# 3. Combine the two to create the population map for your samples; note tab separation is required by STACKS
paste sampleids.txt populations.txt > popmap_geography.txt
```

Second, we will run GSTACKS module in STACKS for alignment and population based summary statistics. Here is the script we will use to do this (run_gstacks.sh):

```bash
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
bam_dir=/n/holyscratch01/hopkins_lab/Chaturvedi/ddrad_sam_april2021/alignment_varcalling/aligned_bamfiles


# Run gstacks to build loci from the aligned paired-end data. We have instructed
# gstacks to remove any PCR duplicates that it finds.
#
gstacks -I $bam_dir -M $bam_dir/popmap_geography.txt --rm-pcr-duplicates -O $stacks -t 8

#
# Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics and 
# smooth the statistics across the genome. Export several output files.
#
populations -P $stacks -M $bam_dir/popmap_geography.txt -r 0.65 --vcf --genepop --fstats --smooth --hwe -t 8
```
