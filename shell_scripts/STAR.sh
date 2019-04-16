#!/bin/bash
#$ -cwd

cd /shelf/apps/User_name/newton/RNAseq/



module load samtools
mkdir star_indicies

# index the genome
/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --runThreadN 16 --limitGenomeGenerateRAM 259760745173 \
--genomeDir ./star_indicies --genomeFastaFiles Gp_Newton_haplotype1.fasta


/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 \
--limitBAMsortRAM 285554136874 --runThreadN 16 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
--outFilterMismatchNmax 7  --outFilterMultimapNmax 5 --outFileNamePrefix Gp \
--readFilesIn  R1.fq.gz R2.fq.gz

