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

###################################################
# Iter 1
ITERNUM=1

#minimap2 -t $THREADS -ax map-pb Gpal_newton_final_strict_scaff_unpolished.fasta ../newton_all_PacBio.fastq.gz > aln.sam
#samtools view -@ $THREADS -S -b -o alnunsorted.bam aln.sam
#wait
#samtools sort -@ $THREADS -o  sorted_mini2.bam alnunsorted.bam
#samtools index sorted_mini2.bam
#pbbamify --input=sorted_mini2.bam --output=sorted_mini.unsorted.pb.bam Gpal_newton_final_strict_scaff_unpolished.fasta /shelf/apps/User_name/newton/mydataset.xml
#samtools sort -@ $THREADS -o sorted_mini.pb.bam sorted_mini.unsorted.pb.bam 
#pbindex  sorted_mini.pb.bam



#/shelf/apps/smrtlink/install/smrtlink-release_6.0.0.47841/bundles/smrttools/install/smrttools-release_6.0.0.47835/private/otherbins/all/bin/arrow  sorted_mini.pb.bam  --diploid --log-file miniarrow.log -j $THREADS --referenceFilename Gpal_newton_final_strict_scaff_unpolished.fasta -o ${PREFIX}.arrow${ITERNUM}.fasta -o ${PREFIX}.arrow${ITERNUM}.gff -o ${PREFIX}.arrow${ITERNUM}.fastq

#rm sorted_mini.unsorted.pb.bam.pbi
#mv sorted_mini.pb.bam sorted_mini.pb1.bam

###################################################
# Iter 2
#minimap2 -t $THREADS -ax map-pb ${PREFIX}.arrow${ITERNUM}.fasta ../newton_all_PacBio.fastq.gz > aln.sam
#samtools view -@ $THREADS -S -b -o alnunsorted.bam aln.sam
#wait 
#samtools sort -@ $THREADS -o  sorted_mini2.bam alnunsorted.bam
#samtools index sorted_mini2.bam
#rm aln.sam alnunsorted.bam
#pbbamify --input=sorted_mini2.bam --output=sorted_mini.unsorted.pb.bam ${PREFIX}.arrow${ITERNUM}.fasta /shelf/apps/User_name/newton/mydataset.xml
#samtools sort -@ $THREADS -o sorted_mini.pb.bam sorted_mini.unsorted.pb.bam 
#pbindex  sorted_mini.pb.bam
#

ITERNUM=2

#/shelf/apps/smrtlink/install/smrtlink-release_6.0.0.47841/bundles/smrttools/install/smrttools-release_6.0.0.47835/private/otherbins/all/bin/arrow  sorted_mini.pb.bam  --diploid --log-file miniarrow.log -j $THREADS --referenceFilename ${PREFIX}.arrow1.fasta -o ${PREFIX}.arrow${ITERNUM}.fasta -o ${PREFIX}.arrow${ITERNUM}.gff -o ${PREFIX}.arrow${ITERNUM}.fastq

#rm sorted_mini.unsorted.pb.bam.pbi
#mv sorted_mini.pb.bam sorted_mini.pb2.bam

###################################################
# Iter 3
#minimap2 -t $THREADS -ax map-pb ${PREFIX}.arrow${ITERNUM}.fasta ../newton_all_PacBio.fastq.gz > aln.sam
#samtools view -@ $THREADS -S -b -o alnunsorted.bam aln.sam
#wait 
#samtools sort -@ $THREADS -o  sorted_mini2.bam alnunsorted.bam
#samtools index sorted_mini2.bam
#rm aln.sam alnunsorted.bam
#pbbamify --input=sorted_mini2.bam --output=sorted_mini.unsorted.pb.bam ${PREFIX}.arrow${ITERNUM}.fasta /shelf/apps/User_name/newton/mydataset.xml
#samtools sort -@ $THREADS -o sorted_mini.pb.bam sorted_mini.unsorted.pb.bam 
#pbindex  sorted_mini.pb.bam

ITERNUM=3

#/shelf/apps/smrtlink/install/smrtlink-release_6.0.0.47841/bundles/smrttools/install/smrttools-release_6.0.0.47835/private/otherbins/all/bin/arrow  sorted_mini.pb.bam  --diploid --log-file miniarrow.log -j $THREADS --referenceFilename ${PREFIX}.arrow2.fasta -o ${PREFIX}.arrow${ITERNUM}.fasta -o ${PREFIX}.arrow${ITERNUM}.gff -o ${PREFIX}.arrow${ITERNUM}.fastq



ITERNUM=2

############################################################################################################################################################
# pilon Illumina polishing-pacbio-assemblies-with-arrow-and-pilon/

##############################################################
##################################################
# pilon 1

#conda activate pilon1.23
THREADS=16
PREFIX=Gp_newton_PB_final

#PILITERNUM=1
#
#minimap2 -t $THREADS -ax map-pb ${PREFIX}.arrow${ITERNUM}.fasta ../newton_all_PacBio.fastq.gz > aln.sam
#samtools view -@ $THREADS -S -b -o alnunsorted.bam aln.sam
#wait
#samtools sort -@ $THREADS -o  sorted_mini2.bam alnunsorted.bam
#samtools index sorted_mini2.bam
#
#bwa index ${PREFIX}.arrow${ITERNUM}.fasta
#bwa mem -t $THREADS ${PREFIX}.arrow${ITERNUM}.fasta /shelf/apps/User_name/newton/DNAseq/R1.fq /shelf/apps/User_name/newton/DNAseq/R2.fq > new.illumina.mapped.sam
#
#samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
#samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam
#rm new.illumina.temp.mapped.bam new.illumina.mapped.sam
#samtools index new.illumina.mapped.bam
#
## this has been through one iteration. Now use the ouptu from iteration 1 for the new error correction. 
## old way = java -Xmx235G -jar ~/scratch/Downloads/pilon-1.22.jar
#pilon --genome ${PREFIX}.arrow${ITERNUM}.fasta --bam new.illumina.mapped.bam --pacbio sorted_mini2.bam --changes --vcf --diploid --threads $THREADS --output Newton_Illumina_pilon_iter1
#
#
###################################################
## pilon 2
#
#PILITERNUM=2
#
#minimap2 -t $THREADS -ax map-pb Newton_Illumina_pilon_iter1.fasta ../newton_all_PacBio.fastq.gz > aln.sam
#samtools view -@ $THREADS -S -b -o alnunsorted.bam aln.sam
#wait
#samtools sort -@ $THREADS -o  sorted_mini2.bam alnunsorted.bam
#samtools index sorted_mini2.bam
#
#bwa index Newton_Illumina_pilon_iter1.fasta
#bwa mem -t $THREADS Newton_Illumina_pilon_iter1.fasta /shelf/apps/User_name/newton/DNAseq/R1.fq /shelf/apps/User_name/newton/DNAseq/R2.fq > new.illumina.mapped.sam
#
#samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
#samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam
#
#samtools index new.illumina.mapped.bam
#rm new.illumina.temp.mapped.bam new.illumina.mapped.sam
#
## this has been through one iteration. Now use the ouptu from iteration 1 for the new error correction. 
#pilon --genome Newton_Illumina_pilon_iter1.fasta --bam new.illumina.mapped.bam --pacbio sorted_mini2.bam --changes --vcf --diploid --threads $THREADS --output Newton_Illumina_pilon_iter2
#
#
###################################################
## pilon 3
#
#PILITERNUM=3
#
#minimap2 -t $THREADS -ax map-pb Newton_Illumina_pilon_iter2.fasta ../newton_all_PacBio.fastq.gz > aln.sam
#samtools view -@ $THREADS -S -b -o alnunsorted.bam aln.sam
#wait
#samtools sort -@ $THREADS -o  sorted_mini2.bam alnunsorted.bam
#samtools index sorted_mini2.bam
#bwa index Newton_Illumina_pilon_iter2.fasta
#bwa mem -t $THREADS Newton_Illumina_pilon_iter2.fasta /shelf/apps/User_name/newton/DNAseq/R1.fq /shelf/apps/User_name/newton/DNAseq/R2.fq > new.illumina.mapped.sam
#
#samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
#samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam
#samtools index new.illumina.mapped.bam
#rm new.illumina.temp.mapped.bam new.illumina.mapped.sam
#
## this has been through one iteration. Now use the ouptu from iteration 1 for the new error correction. 
#pilon --genome Newton_Illumina_pilon_iter2.fasta --bam new.illumina.mapped.bam --pacbio sorted_mini2.bam --changes --vcf --diploid --threads $THREADS --output Newton_Illumina_pilon_iter3
#
#
#
###################################################
#PILITERNUM=4
#
#minimap2 -t $THREADS -ax map-pb Newton_Illumina_pilon_iter3.fasta ../newton_all_PacBio.fastq.gz > aln.sam
#samtools view -@ $THREADS -S -b -o alnunsorted.bam aln.sam
#wait
#samtools sort -@ $THREADS -o  sorted_mini2.bam alnunsorted.bam
#samtools index sorted_mini2.bam
#
#bwa index Newton_Illumina_pilon_iter3.fasta
#bwa mem -t $THREADS Newton_Illumina_pilon_iter3.fasta /shelf/apps/User_name/newton/DNAseq/R1.fq /shelf/apps/User_name/newton/DNAseq/R2.fq > new.illumina.mapped.sam
#
#samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
#samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam
#samtools index new.illumina.mapped.bam
#
#rm new.illumina.temp.mapped.bam new.illumina.mapped.sam
#
## this has been through one iteration. Now use the ouptu from iteration 1 for the new error correction. 
#pilon --genome Newton_Illumina_pilon_iter3.fasta --bam new.illumina.mapped.bam --pacbio sorted_mini2.bam --changes --vcf --diploid --threads $THREADS --output Newton_Illumina_pilon_iter4
#
##################################################
#PILITERNUM=5
#
#minimap2 -t $THREADS -ax map-pb Newton_Illumina_pilon_iter4.fasta ../newton_all_PacBio.fastq.gz > aln.sam
#samtools view -@ $THREADS -S -b -o alnunsorted.bam aln.sam
#wait
#samtools sort -@ $THREADS -o  sorted_mini2.bam alnunsorted.bam
#samtools index sorted_mini2.bam
#
##bwa index Newton_Illumina_pilon_iter4.fasta
#bwa mem -t $THREADS Newton_Illumina_pilon_iter4.fasta /shelf/apps/User_name/newton/DNAseq/R1.fq /shelf/apps/User_name/newton/DNAseq/R2.fq > new.illumina.mapped.sam
#
#samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
#samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam
#
#samtools index new.illumina.mapped.bam
#rm new.illumina.temp.mapped.bam new.illumina.mapped.sam
#
## this has been through one iteration. Now use the ouptu from iteration 1 for the new error correction. 
#pilon --genome Newton_Illumina_pilon_iter4.fasta --bam new.illumina.mapped.bam --pacbio sorted_mini2.bam --changes --vcf --diploid --threads $THREADS --output Newton_Illumina_pilon_iter5
#
###################################################
#PILITERNUM=6
#
#minimap2 -t $THREADS -ax map-pb Newton_Illumina_pilon_iter5.fasta ../newton_all_PacBio.fastq.gz > aln.sam
#samtools view -@ $THREADS -S -b -o alnunsorted.bam aln.sam
#wait
#samtools sort -@ $THREADS -o  sorted_mini2.bam alnunsorted.bam
#samtools index sorted_mini2.bam
#
#bwa index Newton_Illumina_pilon_iter5.fasta
#bwa mem -t $THREADS Newton_Illumina_pilon_iter5.fasta /shelf/apps/User_name/newton/DNAseq/R1.fq /shelf/apps/User_name/newton/DNAseq/R2.fq > new.illumina.mapped.sam
#
#samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
#samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam
#samtools index new.illumina.mapped.bam
#
#rm new.illumina.temp.mapped.bam new.illumina.mapped.sam
#
## this has been through one iteration. Now use the ouptu from iteration 1 for the new error correction. 
#pilon --genome Newton_Illumina_pilon_iter5.fasta --bam new.illumina.mapped.bam --pacbio sorted_mini2.bam --changes --vcf --diploid --threads $THREADS --output Newton_Illumina_pilon_iter6
#
#
#
#################################################
# Haplotype phasing
##################################################
PILITERNUM=7

#conda activate haplo_phase
#
#minimap2 -t $THREADS -ax map-pb Newton_Illumina_pilon_iter5.fasta ../newton_all_PacBio.fastq.gz > aln2.sam
#samtools view -@ $THREADS -S -b -o alnunsorted2.bam aln2.sam
wait
#samtools sort -@ $THREADS -o  sorted_mini3.bam alnunsorted2.bam
#samtools index sorted_mini3.bam

bwa index Newton_Illumina_pilon_iter5.fasta
bwa mem -t $THREADS Newton_Illumina_pilon_iter5.fasta /shelf/apps/User_name/newton/DNAseq/R1_prinseq_good_pm0d.fastq \
/shelf/apps/User_name/newton/DNAseq/R2_prinseq_good_vVZt.fastq > new.illumina.mapped2.sam

samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped2.bam new.illumina.mapped2.sam
samtools sort -@ $THREADS -o  new.illumina.mapped2.bam new.illumina.temp.mapped2.bam
samtools index new.illumina.mapped2.bam


##################################################
# now for haplotype phasing

##################################################
# now for haplotype phasing



#module load freebayes/gitv1_8d2b3a0

freebayes -f Newton_Illumina_pilon_iter5.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.9 -v Newton_illumina.vcf \
--ploidy 2 new.illumina.mapped2.bam

conda activate 
cp sorted_mini3.bam pacbio.bam
whatshap phase --ignore-read-groups -o phased.vcf input.vcf pacbio.bam