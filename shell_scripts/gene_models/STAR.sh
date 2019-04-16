#!/bin/bash
#$ -cwd

#cd /shelf/apps/User_name/newton/RNAseq/

cd /storage/home/users/User_name/newton/gene_models/



module load samtools
#rm -rf star_indicies
#mkdir star_indicies

# index the genome
#~/shelf_apps/apps/STAR-master/bin/Linux_x86_64_static/STAR  --runMode genomeGenerate \
--runThreadN 12 --limitGenomeGenerateRAM 259760745173 \
--genomeDir ./star_indicies --genomeFastaFiles  Gp_Newton_haplotype1.fasta


~/shelf_apps/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ \
--limitGenomeGenerateRAM 405554136874 \
--limitBAMsortRAM 40869285641 \
 --runThreadN 12 --readFilesCommand zcat \
--sjdbGTFfile ./brakerrnaseq/augustus.hints.gtf \
--quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate  \
--outFilterMismatchNmax 7  --outFilterMultimapNmax 5 --outFileNamePrefix Gp \
--readFilesIn  ../RNAseq/R1.fq.gz ../RNAseq/R2.fq.gz


#samtools sort -@ 12  -o sorted.bam Gp*.bam

#samtools index sorted.bam

samtools view -b -f 4 GpAligned.toTranscriptome.out.bam > unmapped.bam

conda activate bedtools

samtools sort -@ 12 -n unmapped.bam aln.qsort


bedtools bamtofastq -i aln.qsort.bam \
                      -fq unmapped_R1.fq \
                      -fq2 unmapped_R2.fq
