#!/bin/bash
#$ -cwd
#set -e

cd /storage/home/users/User_name/newton/final_genome2/RNAseq_mapping/

conda activate python27

# HTseq counts

#/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 16 \
 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 5 \
 --outFileNamePrefix Gp_MALE_rep2 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_MALE_rep2_ERR202435_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_MALE_rep2_ERR202435_paired_2.fastq.gz 
#
#
#samtools sort -@ 16 -o Gp_MALE_rep2*.bam Gp_MALE_rep2Aligned.sortedByCoord.out.bam
#samtools index Gp_MALE_rep2*.bam 
 
cd /storage/home/users/User_name/newton/final_genome2/RNAseq_mapping/

htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_MALE_rep2Aligned.sortedByCoord.out.bam Gpal_newton_newton.gff3 > Gp_MALE_rep2_genes.counts

echo "	Gp_7DPI_rep1	Gp_7DPI_rep2	Gp_14DPI_rep1	Gp_14DPI_rep2	Gp_21DPI_rep1	Gp_21DPI_rep2	Gp_28DPI_rep1	Gp_28DPI_rep2	Gp_35DPI_rep1	Gp_35DPI_rep2	Gp_EGG_rep1	Gp_EGG_rep2	Gp_J2_rep1	Gp_J2_rep2	Gp_J2_rep3	Gp_MALE_rep1	Gp_MALE_rep2" > Gp_genes.counts.matrix
# merge these into a file
FILES=$(ls -t -v *.counts | tr '\n' ' ')
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES | grep -v "__" >> Gp_genes.counts.matrix
