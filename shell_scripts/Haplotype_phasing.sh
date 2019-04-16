#!/bin/bash
#$ -cwd

cd /shelf/apps/User_name/newton

haplotype phasing

https://whatshap.readthedocs.io/en/latest/guide.html#recommended-workflow

https://whatshap.readthedocs.io/en/latest/guide.html#subcommands

User_name@marvin:~ > conda activate phasing
(phasing) User_name@marvin:~ > whatshap -h
usage: whatshap [-h] [--version] [--debug]
                {phase,stats,compare,hapcut2vcf,unphase,haplotag,genotype} ...

positional arguments:
  {phase,stats,compare,hapcut2vcf,unphase,haplotag,genotype}
    phase               Phase variants in a VCF with the WhatsHap algorithm
    stats               Print phasing statistics of a single VCF file
    compare             Compare two or more phasings
    hapcut2vcf          Convert hapCUT output format to VCF
    unphase             Remove phasing information from a VCF file
    haplotag            Tag reads by haplotype
    genotype            Genotype variants

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --debug               Print debug messages
(phasing) User_name@marvin:~ >

#conda activate pilon1.23
THREADS=16
PREFIX=Gp_newton_PB_final


##################################################
PILITERNUM=6

bwa index Newton_Illumina_pilon_iter6.fasta
bwa mem -t $THREADS Newton_Illumina_pilon_iter6.fasta /shelf/apps/User_name/newton/DNAseq/R1.fq /shelf/apps/User_name/newton/DNAseq/R2.fq > new.illumina.mapped.sam

samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam
samtools index new.illumina.mapped.bam

minimap2 -t $THREADS -ax map-pb Newton_Illumina_pilon_iter6.fasta ../newton_all_PacBio.fastq.gz > aln.sam
samtools view -@ $THREADS -S -b -o alnunsorted.bam aln.sam
wait
samtools sort -@ $THREADS -o  sorted_mini2.bam alnunsorted.bam
samtools index sorted_mini2.bam

module load freebayes/gitv1_8d2b3a0

freebayes -f Newton_Illumina_pilon_iter6.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.9 -v Newton_illumina.vcf --ploidy 2 new.illumina.mapped.bam

cp sorted_mini2.bam pacbio.bam
whatshap phase -o phased.vcf input.vcf pacbio.bam
--ignore-read-groups


bgzip phased.vcf
tabix phased.vcf.gz
bcftools consensus -H 1 -f reference.fasta phased.vcf.gz > haplotype1.fasta
bcftools consensus -H 2 -f reference.fasta phased.vcf.gz > haplotype2.fasta

whatshap stats --gtf=phased.gtf phased.vcf

