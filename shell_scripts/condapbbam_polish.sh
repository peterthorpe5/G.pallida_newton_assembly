#!/bin/bash
#$ -cwd
set -e

cd /shelf/apps/User_name/newton

#11) POLISHING. Not yet achived. The default approach (a different mapper will be chosen due to the slow run time of BLASR)
#https:/flowersoftheocean.wordpress.com/2018/04/16/polishing-pacbio-assemblies-with-arrow-and-pilon/
# blasr , although slow has some good option, such as seed. Alignment lenght etc ... see below. This takes weeks to run. 
# Polishing needs to be repeated 3 or more time, this adding months onto assembly-process time.

conda activate pb_toolkit


THREADS=16
PREFIX=Gp_newton_conda_PB



#singularity run --app arrow /shelf/apps/smrt-tools.img sorted_min2i.pb.bam  --diploid --log-file miniarrow.log -j $THREADS --referenceFilename ${PREFIX}.arrow1.fasta -o ${PREFIX}.arrow${ITERNUM}.fasta -o ${PREFIX}.arrow${ITERNUM}.gff -o ${PREFIX}.arrow${ITERNUM}.fastq
#singularity run --app arrow /shelf/apps/smrt-tools.img  sorted_min2i.pb.bam  --diploid --log-file miniarrow.log -j 1 --referenceFilename Gpal_newton_final_unpolished.fasta -o ${PREFIX}.arrow${ITERNUM}.fasta -o ${PREFIX}.arrow${ITERNUM}.gff -o ${PREFIX}.arrow${ITERNUM}.fastq



###################################################
# Iter 1
ITERNUM=1

##minimap2 -t $THREADS -ax map-pb Gpal_newton_final_unpolished.fasta newton_all_PacBio.fastq.gz > aln2.sam
#samtools view -@ $THREADS -S -b -o alnunsorted2.bam aln2.sam
wait
#samtools sort -@ $THREADS -o  sorted_mini3.bam alnunsorted2.bam
#samtools index sorted_mini3.bam
#pbbamify --input=sorted_mini3.bam --output=sorted_mini.unsorted2.pb.bam Gpal_newton_final_unpolished.fasta /shelf/apps/User_name/newton/mydataset.xml
#samtools sort -@ $THREADS -o sorted_min2i.pb.bam sorted_mini.unsorted2.pb.bam 
#pbindex  sorted_min2i.pb.bam


arrow sorted_mini.pb.bam  --diploid --log-file miniarrow.log -j $THREADS --referenceFilename Gpal_newton_final_unpolished.fasta -o ${PREFIX}.arrow${ITERNUM}.fasta -o ${PREFIX}.arrow${ITERNUM}.gff -o ${PREFIX}.arrow${ITERNUM}.fastq


rm sorted_mini.unsorted2.pb.bam.pbi
mv sorted_min2i.pb.bam sorted_mini.pb1.bam

###################################################
# Iter 2
minimap2 -t $THREADS -ax map-pb ${PREFIX}.arrow${ITERNUM}.fasta newton_all_PacBio.fastq.gz > aln2.sam
samtools view -@ $THREADS -S -b -o alnunsorted2.bam aln2.sam
wait 
samtools sort -@ $THREADS -o  sorted_mini3.bam alnunsorted2.bam
samtools index sorted_mini3.bam
rm aln2.sam alnunsorted2.bam
pbbamify --input=sorted_mini3.bam --output=sorted_mini.unsorted2.pb.bam ${PREFIX}.arrow${ITERNUM}.fasta /shelf/apps/User_name/newton/mydataset.xml
samtools sort -@ $THREADS -o sorted_min2i.pb.bam sorted_mini.unsorted2.pb.bam 
pbindex  sorted_min2i.pb.bam


ITERNUM=2
arrow  sorted_min2i.pb.bam  --diploid --log-file miniarrow.log -j $THREADS --referenceFilename ${PREFIX}.arrow1.fasta -o ${PREFIX}.arrow${ITERNUM}.fasta -o ${PREFIX}.arrow${ITERNUM}.gff -o ${PREFIX}.arrow${ITERNUM}.fastq

rm sorted_mini.unsorted2.pb.bam.pbi
mv sorted_min2i.pb.bam sorted_mini.pb2.bam

###################################################
# Iter 3
minimap2 -t $THREADS -ax map-pb ${PREFIX}.arrow${ITERNUM}.fasta newton_all_PacBio.fastq.gz > aln2.sam
samtools view -@ $THREADS -S -b -o alnunsorted2.bam aln2.sam
wait 
samtools sort -@ $THREADS -o  sorted_mini3.bam alnunsorted2.bam
samtools index sorted_mini3.bam
rm aln2.sam alnunsorted2.bam
pbbamify --input=sorted_mini3.bam --output=sorted_mini.unsorted2.pb.bam ${PREFIX}.arrow${ITERNUM}.fasta /shelf/apps/User_name/newton/mydataset.xml
samtools sort -@ $THREADS -o sorted_min2i.pb.bam sorted_mini.unsorted2.pb.bam 
pbindex  sorted_min2i.pb.bam

ITERNUM=3
arrow  sorted_min2i.pb.bam  --diploid --log-file miniarrow.log -j $THREADS --referenceFilename ${PREFIX}.arrow2.fasta -o ${PREFIX}.arrow${ITERNUM}.fasta -o ${PREFIX}.arrow${ITERNUM}.gff -o ${PREFIX}.arrow${ITERNUM}.fastq





############################################################################################################################################################
# pilon Illumina polishing-pacbio-assemblies-with-arrow-and-pilon/

##############################################################
##################################################
# pilon 1

PILITERNUM=1

bwa index ${PREFIX}.arrow${ITERNUM}.fasta
bwa mem -t $THREADS ${PREFIX}.arrow${ITERNUM}.fasta /shelf/apps/User_name/newton/DNAseq/R1.fq /shelf/apps/User_name/newton/DNAseq/R2.fq > new.illumina.mapped.sam

samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam
rm new.illumina.temp.mapped.bam new.illumina.mapped.sam
samtools index new.illumina.mapped.bam

# this has been through one iteration. Now use the ouptu from iteration 1 for the new error correction. 
java -Xmx235G -jar ~/scratch/Downloads/pilon-1.22.jar --genome ${PREFIX}.arrow${ITERNUM}.fasta --bam new.illumina.mapped.bam --changes --vcf --diploid --threads $THREADS --output Newton_Illumina_pilon_iter1


##################################################
# pilon 2

PILITERNUM=2

bwa index Newton_Illumina_pilon_iter1.fasta
bwa mem -t $THREADS Newton_Illumina_pilon_iter1.fasta /shelf/apps/User_name/newton/DNAseq/R1.fq /shelf/apps/User_name/newton/DNAseq/R2.fq > new.illumina.mapped.sam

samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam

samtools index new.illumina.mapped.bam
rm new.illumina.temp.mapped.bam new.illumina.mapped.sam

# this has been through one iteration. Now use the ouptu from iteration 1 for the new error correction. 
java -Xmx220G -jar ~/scratch/Downloads/pilon-1.22.jar --genome Newton_Illumina_pilon_iter1.fasta --bam new.illumina.mapped.bam --changes --vcf --diploid --threads $THREADS --output Newton_Illumina_pilon_iter2


##################################################
# pilon 3

PILITERNUM=3

bwa index Newton_Illumina_pilon_iter2.fasta
bwa mem -t $THREADS Newton_Illumina_pilon_iter2.fasta /shelf/apps/User_name/newton/DNAseq/R1.fq /shelf/apps/User_name/newton/DNAseq/R2.fq > new.illumina.mapped.sam

samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam
samtools index new.illumina.mapped.bam
rm new.illumina.temp.mapped.bam new.illumina.mapped.sam

# this has been through one iteration. Now use the ouptu from iteration 1 for the new error correction. 
java -Xmx220G -jar ~/scratch/Downloads/pilon-1.22.jar --genome Newton_Illumina_pilon_iter2.fasta --bam new.illumina.mapped.bam --changes --vcf --diploid --threads $THREADS --output Newton_Illumina_pilon_iter3



##################################################
PILITERNUM=4

bwa index Newton_Illumina_pilon_iter3.fasta
bwa mem -t $THREADS Newton_Illumina_pilon_iter3.fasta /shelf/apps/User_name/newton/DNAseq/R1.fq /shelf/apps/User_name/newton/DNAseq/R2.fq > new.illumina.mapped.sam

samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam
samtools index new.illumina.mapped.bam

rm new.illumina.temp.mapped.bam new.illumina.mapped.sam

# this has been through one iteration. Now use the ouptu from iteration 1 for the new error correction. 
java -Xmx230G -jar ~/scratch/Downloads/pilon-1.22.jar --genome Newton_Illumina_pilon_iter3.fasta --bam new.illumina.mapped.bam --changes --vcf --diploid --threads $THREADS --output Newton_Illumina_pilon_iter4

##################################################
PILITERNUM=5

bwa index Newton_Illumina_pilon_iter4.fasta
bwa mem -t $THREADS Newton_Illumina_pilon_iter4.fasta /shelf/apps/User_name/newton/DNAseq/R1.fq /shelf/apps/User_name/newton/DNAseq/R2.fq > new.illumina.mapped.sam

samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam

samtools index new.illumina.mapped.bam
rm new.illumina.temp.mapped.bam new.illumina.mapped.sam

# this has been through one iteration. Now use the ouptu from iteration 1 for the new error correction. 
java -Xmx230G -jar ~/scratch/Downloads/pilon-1.22.jar --genome Newton_Illumina_pilon_iter4.fasta --bam new.illumina.mapped.bam --changes --vcf --diploid --threads $THREADS --output Newton_Illumina_pilon_iter5

##################################################
PILITERNUM=6

bwa index Newton_Illumina_pilon_iter5.fasta
bwa mem -t $THREADS Newton_Illumina_pilon_iter5.fasta /shelf/apps/User_name/newton/DNAseq/R1.fq /shelf/apps/User_name/newton/DNAseq/R2.fq > new.illumina.mapped.sam

samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped.bam new.illumina.mapped.sam
samtools sort -@ $THREADS -o  new.illumina.mapped.bam new.illumina.temp.mapped.bam
samtools index new.illumina.mapped.bam

rm new.illumina.temp.mapped.bam new.illumina.mapped.sam

# this has been through one iteration. Now use the ouptu from iteration 1 for the new error correction. 
java -Xmx230G -jar ~/scratch/Downloads/pilon-1.22.jar --genome Newton_Illumina_pilon_iter5.fasta --bam new.illumina.mapped.bam --changes --vcf --diploid --threads $THREADS --output Newton_Illumina_pilon_iter6

