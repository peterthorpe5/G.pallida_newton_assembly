#!/bin/bash
#$ -cwd
#set -e


cd $HOME/newton/final_genome2

conda activate python27

#module load samtools
rm -rf star_indicies
mkdir star_indicies

# index the genome
/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --runThreadN 32  --limitGenomeGenerateRAM 259760745173 \
--genomeDir ./star_indicies \
--genomeFastaFiles Gp_Newton_haplotype1.fasta

########################################################################
 #star mapping
 #star mapping

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_14DPI_rep1 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_14DPI_rep1_ERR202423_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_14DPI_rep1_ERR202423_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_14DPI_rep1*.bam Gp_14DPI_rep1Aligned.sortedByCoord.out.bam
samtools index Gp_14DPI_rep1*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_14DPI_rep2 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_14DPI_rep2_ERR407001_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_14DPI_rep2_ERR407001_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_14DPI_rep2*.bam Gp_14DPI_rep2Aligned.sortedByCoord.out.bam
samtools index Gp_14DPI_rep2*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_21DPI_rep1 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_21DPI_rep1_ERR406997_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_21DPI_rep1_ERR406997_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_21DPI_rep1*.bam Gp_21DPI_rep1Aligned.sortedByCoord.out.bam
samtools index Gp_21DPI_rep1*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_21DPI_rep2 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_21DPI_rep2_ERR202426_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_21DPI_rep2_ERR202426_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_21DPI_rep2*.bam Gp_21DPI_rep2Aligned.sortedByCoord.out.bam
samtools index Gp_21DPI_rep2*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_28DPI_rep1 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_28DPI_rep1_ERR406995_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_28DPI_rep1_ERR406995_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_28DPI_rep1*.bam Gp_28DPI_rep1Aligned.sortedByCoord.out.bam
samtools index Gp_28DPI_rep1*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_28DPI_rep2 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_28DPI_rep2_ERR202427_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_28DPI_rep2_ERR202427_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_28DPI_rep2*.bam Gp_28DPI_rep2Aligned.sortedByCoord.out.bam
samtools index Gp_28DPI_rep2*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_35DPI_rep1 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_35DPI_rep1_ERR406998_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_35DPI_rep1_ERR406998_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_35DPI_rep1*.bam Gp_35DPI_rep1Aligned.sortedByCoord.out.bam
samtools index Gp_35DPI_rep1*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_35DPI_rep2 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_35DPI_rep2_ERR202428_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_35DPI_rep2_ERR202428_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_35DPI_rep2*.bam Gp_35DPI_rep2Aligned.sortedByCoord.out.bam
samtools index Gp_35DPI_rep2*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_7DPI_rep1 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_7DPI_rep1_ERR202425_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_7DPI_rep1_ERR202425_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_7DPI_rep1*.bam Gp_7DPI_rep1Aligned.sortedByCoord.out.bam
samtools index Gp_7DPI_rep1*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_7DPI_rep2 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_7DPI_rep2_ERR202436_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_7DPI_rep2_ERR202436_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_7DPI_rep2*.bam Gp_7DPI_rep2Aligned.sortedByCoord.out.bam
samtools index Gp_7DPI_rep2*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_EGG_rep1 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_EGG_rep1_ERR407000_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_EGG_rep1_ERR407000_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_EGG_rep1*.bam Gp_EGG_rep1Aligned.sortedByCoord.out.bam
samtools index Gp_EGG_rep1*.bam

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_EGG_rep2 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_EGG_rep2_ERR202430_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_EGG_rep2_ERR202430_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_EGG_rep2*.bam Gp_EGG_rep2Aligned.sortedByCoord.out.bam
samtools index Gp_EGG_rep2*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_J2_rep1 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_J2_rep1_ERR406996_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_J2_rep1_ERR406996_paired_2.fastq.gz 

 # 14
#samtools sort -@ 8 -o Gp_J2_rep1*.bam Gp_J2_rep1Aligned.sortedByCoord.out.bam
samtools index Gp_J2_rep1*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 12 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_J2_rep2 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_J2_rep2_ERR202424_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_J2_rep2_ERR202424_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_J2_rep2*.bam Gp_J2_rep2Aligned.sortedByCoord.out.bam
samtools index Gp_J2_rep2*.bam 



/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 12 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_J2_rep3 --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_J2_rep3_ERR202429_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_J2_rep3_ERR202429_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_J2_rep3*.bam Gp_J2_rep3Aligned.sortedByCoord.out.bam
samtools index Gp_J2_rep3*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 6 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 --outFileNamePrefix Gp_MALE_rep1_c --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_MALE_rep1_ERR202422_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_MALE_rep1_ERR202422_paired_2.fastq.gz 


#samtools sort -@ 8 -o Gp_MALE_rep1*c.bam Gp_MALE_rep1Aligned.sortedByCoord.out.bam
samtools index Gp_MALE_rep1*.bam 

/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 32  \
 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 3 \
 --outFileNamePrefix Gp_MALE_rep2_b --readFilesIn /storage/home/users/User_name/project/Gp_lifestages/Gp_MALE_rep2_ERR202435_paired_1.fastq.gz /storage/home/users/User_name/project/Gp_lifestages/Gp_MALE_rep2_ERR202435_paired_2.fastq.gz 
#
#
##samtools sort -@ 8 -o Gp_MALE_rep2*.bam Gp_MALE_rep2Aligned.sortedByCoord.out.bam
#samtools index Gp_MALE_rep2*.bam 

####################################################################################################
# HTseq counts

#conda activate python27
#htseq-count --quiet -f bam --type gene --idattr ID Gp_14DPI_rep1*.bam ../augustus.hints.gff3 > Gp_14DPI_rep1_genes.counts
#htseq-count --quiet -f bam --type gene --idattr ID Gp_14DPI_rep2*.bam ../augustus.hints.gff3 > Gp_14DPI_rep2_genes.counts
#htseq-count --quiet -f bam --type gene --idattr ID Gp_21DPI_rep1*.bam ../augustus.hints.gff3 > Gp_21DPI_rep1_genes.counts


#htseq-count --quiet -f bam --type gene --idattr ID Gp_21DPI_rep2*.bam ../augustus.hints.gff3 > Gp_21DPI_rep2_genes.counts
#htseq-count --quiet -f bam --type gene --idattr ID Gp_28DPI_rep1*.bam ../augustus.hints.gff3 > Gp_28DPI_rep1_genes.counts

#htseq-count --quiet -f bam --type gene --idattr ID Gp_28DPI_rep2*.bam ../augustus.hints.gff3 > Gp_28DPI_rep2_genes.counts
#htseq-count --quiet -f bam --type gene --idattr ID Gp_35DPI_rep1*.bam ../augustus.hints.gff3 > Gp_35DPI_rep1_genes.counts

#htseq-count --quiet -f bam --type gene --idattr ID Gp_35DPI_rep2*.bam ../augustus.hints.gff3 > Gp_35DPI_rep2_genes.counts
#htseq-count --quiet -f bam --type gene --idattr ID Gp_7DPI_rep1*.bam ../augustus.hints.gff3 > Gp_7DPI_rep1_genes.counts

#htseq-count --quiet -f bam --type gene --idattr ID Gp_7DPI_rep2*.bam ../augustus.hints.gff3 > Gp_7DPI_rep2_genes.counts
#htseq-count --quiet -f bam --type gene --idattr ID Gp_EGG_rep1*.bam ../augustus.hints.gff3 > Gp_EGG_rep1_2_genes.counts

#htseq-count --quiet -f bam --type gene --idattr ID Gp_EGG_rep2*.bam ../augustus.hints.gff3 > Gp_EGG_rep2_genes.counts
#htseq-count --quiet -f bam --type gene --idattr ID Gp_J2_rep1*.bam ../augustus.hints.gff3 > Gp_J2_rep1_genes.counts
#htseq-count --quiet -f bam --type gene --idattr ID Gp_J2_rep2*.bam ../augustus.hints.gff3 > Gp_J2_rep2_genes.counts

#htseq-count --quiet -f bam --type gene --idattr ID Gp_J2_rep3*.bam ../augustus.hints.gff3 > Gp_J2_rep3_genes.counts
#htseq-count --quiet -f bam --type gene --idattr ID Gp_MALE_rep1*.bam ../augustus.hints.gff3 > Gp_MALE_rep1_genes.counts
htseq-count --quiet -f bam --type gene --idattr ID Gp_MALE_rep2*.bam ../augustus.hints.gff3 > Gp_MALE_rep2_1_genes.counts

echo "\tGp_7DPI_rep1\tGp_7DPI_rep2\tGp_14DPI_rep1\tGp_14DPI_rep2\tGp_21DPI_rep1\tGp_21DPI_rep2\tGp_28DPI_rep1\tGp_28DPI_rep2\tGp_35DPI_rep1\tGp_35DPI_rep2\tGp_EGG_rep1\tGp_EGG_rep2\tGp_J2_rep1\tGp_J2_rep2\tGp_J2_rep3\tGp_MALE_rep1\tGp_MALE_rep2" > Gp_genes.counts.matrix
# merge these into a file
FILES=$(ls -t -v *.counts | tr '\n' ' ')
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES | grep -v "__" >> Gp_genes.counts.matrix
