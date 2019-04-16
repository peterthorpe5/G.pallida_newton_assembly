#!/bin/bash
#$ -cwd
set -e

cd /shelf/apps/User_name/newton/stricter_scaff_polish


#11) POLISHING. Not yet achived. The default approach (a different mapper will be chosen due to the slow run time of BLASR)
#https:/flowersoftheocean.wordpress.com/2018/04/16/polishing-pacbio-assemblies-with-arrow-and-pilon/
# blasr , although slow has some good option, such as seed. Alignment lenght etc ... see below. This takes weeks to run. 
# Polishing needs to be repeated 3 or more time, this adding months onto assembly-process time.

#conda activate pb_toolkit

THREADS=16
PREFIX=Gp_newton_PB_final


PILITERNUM=7

#conda activate haplo_phase
#
#minimap2 -t $THREADS -ax map-pb Newton_Illumina_pilon_iter5.fasta ../newton_all_PacBio.fastq.gz > aln2.sam
#samtools view -@ $THREADS -S -b -o alnunsorted2.bam aln2.sam
wait
#samtools sort -@ $THREADS -o  sorted_mini3.bam alnunsorted2.bam
#samtools index sorted_mini3.bam

#bwa index Gp_Newton.fasta
#bwa mem -t $THREADS Gp_Newton.fasta /shelf/apps/User_name/newton/DNAseq/R1_prinseq_good_pm0d.fastq \
#/shelf/apps/User_name/newton/DNAseq/R2_prinseq_good_vVZt.fastq > new.illumina.mapped2.sam

#samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped2.bam new.illumina.mapped2.sam
#samtools sort -@ $THREADS -o  new.illumina.mapped2.bam new.illumina.temp.mapped2.bam
#samtools index new.illumina.mapped2.bam


##################################################
# now for haplotype phasing

##################################################
# now for haplotype phasing



#module load freebayes/gitv1_8d2b3a0

#freebayes -f Gp_Newton.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.9 -v Newton_illumina.vcf \
--ploidy 2 new.illumina.mapped2.bam

conda activate haplo_phase
cp sorted_mini3.bam pacbio.bam
samtools index pacbio.bam
whatshap phase --ignore-read-groups -o phased.vcf Newton_illumina.vcf pacbio.bam




