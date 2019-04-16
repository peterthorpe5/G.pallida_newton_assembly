#!/bin/bash
#$ -cwd

cd /shelf/apps/User_name/newton/RNAseq/

sleep 43200

mkdir arrow_indicies3

# index the genome
/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --runThreadN 12 --limitGenomeGenerateRAM 259760745173 \
--genomeDir ./arrow_indicies3 --genomeFastaFiles Gpal_newton_final_strict_scaff_unpolished.fasta 


# for all haplotypes!! allow more mapping.
/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR  --outSAMunmapped None --outReadsUnmapped Unmapped.out.mate --genomeDir arrow_indicies3/ \
 --runThreadN 12 --readFilesCommand zcat --outFilterMismatchNmax 10  \
 --outFilterMultimapNmax 5 --outFileNamePrefix FINAL --readFilesIn  R1.fq.gz R2.fq.gz
 
##rm -rf arrow_indicies2
##mkdir arrow_indicies2
##
##
### index the genome
##shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --runThreadN 4 --limitGenomeGenerateRAM 259760745173 \
##-genomeDir ./arrow_indicies2 --genomeFastaFiles scaffolds.fasta
##
##
## for all haplotypes!! allow more mapping.
##shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR  --outSAMunmapped None --outReadsUnmapped Unmapped.out.mate --genomeDir arrow_indicies2/ \
##--runThreadN 8 --readFilesCommand zcat --outFilterMismatchNmax 10  \
##--outFilterMultimapNmax 5 --outFileNamePrefix scaff_V3 --readFilesIn  R1.fq.gz R2.fq.gz
##
##pigz -p 16 *.fq
##
##
########################################################################################################
##samtools sort -@ 4 -o GpAligned.sortedByCoord.out.bam_sorted.bam GpAligned.sortedByCoord.out.bam
##samtools sort -@ 4 GpAligned.sortedByCoord.out.bam GpAligned.sortedByCoord.out.bam_sorted 
##
### index the bam files
##samtools index GpAligned.sortedByCoord.out.bam_sorted.bam 
##
##samtools index GpAligned.sortedByCoord.out.bam
##
##bam2hints --intronsonly --in=GpAligned.sortedByCoord.out.bam --out=hints.gff
## merge these into a file
###!/bin/bash
###$ -cwd
###-l hostname="n09-08-144-biggus"
##cd /home/pt40963/scratch/nematode/G.pallida/RNAseq
##
###bowtie2-build -f Phytophthora_infestans.ASM14294v1.31.fa Pi
###cp /mnt/shared/projects/nematodes/20180501_Heterodera_schachtii/RNASEQ/* ./RNAseq/
##################################
### T30-4
###cat /mnt/shared/projects/nematodes/Globodera_pallida/ERP001236_transcriptome_data/Gp_lifestages/*_1.fastq.gz > r1.fq.gz
##
###cat /mnt/shared/projects/nematodes/Globodera_pallida/ERP001236_transcriptome_data/Gp_lifestages/*_2.fastq.gz > r2.fq.gz
##
###cat /mnt/shared/projects/nematodes/Globodera_pallida/ERP001236_transcriptome_data/Gp_populations/Gp_Newton_1997_J2_ERR202432_1.fastq.gz R1.fq.gz > r1.fq.gz
###cat /mnt/shared/projects/nematodes/Globodera_pallida/ERP001236_transcriptome_data/Gp_populations/Gp_Newton_1997_J2_ERR202432_2.fastq.gz R2.fq.gz > r2.fq.gz
##
###java -jar /home/pt40963/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 4 -phred33 r1.fq.gz r2.fq.gz R1.fq.gz crap_R1_unpaired.fq.gz R2.fq.gz crap_R2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:67 
##
###rm crap_*
##
###fix fastq with bbmap
##
###repair.sh in1=R1.fq.gz in2=R2.fq.gz out1=fixed1.fq.gz out2=fixed2.fq.gz outsingle=singletons.fq.gz
## 
###cp /mnt/shared/projects/nematodes/Globodera_pallida/20180507_Gp_Newton/assembly/canu16_ALL_SMRT/GBN_All4.contigs.fasta ./
##
###mkdir star_indicies
##
### index the genome
##~/scratch/Downloads/STAR/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --runThreadN 8 --limitGenomeGenerateRAM 259760745173 --genomeDir ./star_indicies --genomeFastaFiles Newton_Illumina_pilon_iter6.fasta
##
### for all haplotypes!! allow more mapping.
##~/scratch/Downloads/STAR/bin/Linux_x86_64_static/STAR --outReadsUnmapped Unmapped.out.mate1 Unmapped.out.mate2 --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 5 --outFileNamePrefix Gp --readFilesIn  R1.fq.gz R2.fq.gz
##
##
###samtools sort -@ 4 -o GpAligned.sortedByCoord.out.bam_sorted.bam GpAligned.sortedByCoord.out.bam
###samtools sort -@ 4 GpAligned.sortedByCoord.out.bam GpAligned.sortedByCoord.out.bam_sorted 
##
#### index the bam files
##samtools index GpAligned.sortedByCoord.out.bam_sorted.bam 
##
##samtools index GpAligned.sortedByCoord.out.bam
##
##bam2hints --intronsonly --in=GpAligned.sortedByCoord.out.bam --out=hints.gff
## merge these into a file