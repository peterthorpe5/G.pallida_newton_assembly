#!/bin/bash
#$ -cwd
cd /storage/home/users/User_name/project/GPAL_Lindley
THREADS=8

java -jar /shelf/training/Trimmomatic-0.38/trimmomatic-0.38.jar PE -summary trim_summary.txt \
-threads $THREADS -phred33 Gp_Lindley_ERR114517_1.fastq.gz Gp_Lindley_ERR114517_2.fastq.gz \
R1_paired.fastq.gz \
R1_unpaired.fastq.gz \
R2_paired.fastq.gz \
R2_unpaired.fastq.gz \
ILLUMINACLIP:/shelf/training/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:45 

cd /storage/home/users/User_name/Lindley

#haplotype phasing


bwa index Gp_Newton_haplotype1.fasta
bwa mem -t $THREADS Gp_Newton_haplotype1.fasta /storage/home/users/User_name/project/GPAL_Lindley/R1_paired.fastq.gz \
 /storage/home/users/User_name/project/GPAL_Lindley/R2_paired.fastq.gz > new.illumina.mapped.sam

samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam
samtools index new.illumina.mapped.bam


module load freebayes/gitv1_8d2b3a0

freebayes -f Gp_Newton_haplotype1.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction \
0.9 -v Newton_illumina.vcf --ploidy 2 new.illumina.mapped.bam

#cp sorted_mini2.bam pacbio.bam
#whatshap phase -o phased.vcf input.vcf pacbio.bam
#--ignore-read-groups


bgzip phased.vcf
tabix phased.vcf.gz
bcftools consensus -H 1 -f reference.fasta phased.vcf.gz > haplotype1.fasta
bcftools consensus -H 2 -f reference.fasta phased.vcf.gz > haplotype2.fasta

whatshap stats --gtf=phased.gtf phased.vcf

