#!/bin/bash
#$ -cwd

cd /shelf/apps/User_name/newton/RNAseq/



#mkdir star_indicies

# index the genome
#/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --runThreadN 8 --limitGenomeGenerateRAM 259760745173 --genomeDir ./star_indicies \
#--genomeFastaFiles scaffolds_gapfilled_FINAL.fasta

# for all haplotypes!! allow more mapping.
#/shelf/apps/User_name/apps/STAR-master/bin/Linux_x86_64_static/STAR --outReadsUnmapped fastq --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 \
#--limitBAMsortRAM 285554136874 --runThreadN 12 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 5 \
#--outFileNamePrefix Gpbraker --readFilesIn  R1.fq.gz R2.fq.gz


########


#######################################################################################################
#samtools sort -@ 4 -o GpAligned.sortedByCoord.out.bam_sorted.bam GpAligned.sortedByCoord.out.bam
#samtools sort -@ 4 GpAligned.sortedByCoord.out.bam GpAligned.sortedByCoord.out.bam_sorted 

## index the bam files
#samtools index GpbrakerAligned.sortedByCoord.out.bam

#bam2hints --intronsonly --in=GpbrakerAligned.sortedByCoord.out.bam --out=hints.gff
# merge these into a file
##!/bin/bash
##$ -cwd
##-l hostname="n09-08-144-biggus"
#cd /home/pt40963/scratch/nematode/G.pallida/RNAseq
#
##bowtie2-build -f Phytophthora_infestans.ASM14294v1.31.fa Pi
##cp /mnt/shared/projects/nematodes/20180501_Heterodera_schachtii/RNASEQ/* ./RNAseq/
#################################
## T30-4
##cat /mnt/shared/projects/nematodes/Globodera_pallida/ERP001236_transcriptome_data/Gp_lifestages/*_1.fastq.gz > r1.fq.gz
#
##cat /mnt/shared/projects/nematodes/Globodera_pallida/ERP001236_transcriptome_data/Gp_lifestages/*_2.fastq.gz > r2.fq.gz
#
##cat /mnt/shared/projects/nematodes/Globodera_pallida/ERP001236_transcriptome_data/Gp_populations/Gp_Newton_1997_J2_ERR202432_1.fastq.gz R1.fq.gz > r1.fq.gz
##cat /mnt/shared/projects/nematodes/Globodera_pallida/ERP001236_transcriptome_data/Gp_populations/Gp_Newton_1997_J2_ERR202432_2.fastq.gz R2.fq.gz > r2.fq.gz
#
##java -jar /home/pt40963/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 4 -phred33 r1.fq.gz r2.fq.gz R1.fq.gz crap_R1_unpaired.fq.gz R2.fq.gz crap_R2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:67 
#
##rm crap_*
#
##fix fastq with bbmap
#
##repair.sh in1=R1.fq.gz in2=R2.fq.gz out1=fixed1.fq.gz out2=fixed2.fq.gz outsingle=singletons.fq.gz
# 
##cp /mnt/shared/projects/nematodes/Globodera_pallida/20180507_Gp_Newton/assembly/canu16_ALL_SMRT/GBN_All4.contigs.fasta ./
#
##mkdir star_indicies
#
## index the genome
#~/scratch/Downloads/STAR/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --runThreadN 8 --limitGenomeGenerateRAM 259760745173 --genomeDir ./star_indicies --genomeFastaFiles Newton_Illumina_pilon_iter6.fasta
#
## for all haplotypes!! allow more mapping.
#~/scratch/Downloads/STAR/bin/Linux_x86_64_static/STAR --outReadsUnmapped Unmapped.out.mate1 Unmapped.out.mate2 --genomeDir star_indicies/ --limitGenomeGenerateRAM 285554136874 --limitBAMsortRAM 285554136874 --runThreadN 8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 5 --outFileNamePrefix Gp --readFilesIn  R1.fq.gz R2.fq.gz
#
#
##samtools sort -@ 4 -o GpAligned.sortedByCoord.out.bam_sorted.bam GpAligned.sortedByCoord.out.bam
##samtools sort -@ 4 GpAligned.sortedByCoord.out.bam GpAligned.sortedByCoord.out.bam_sorted 
#
### index the bam files
#samtools index GpAligned.sortedByCoord.out.bam_sorted.bam 
#
#samtools index GpAligned.sortedByCoord.out.bam
#
#bam2hints --intronsonly --in=GpAligned.sortedByCoord.out.bam --out=hints.gff
# merge these into a file


#################
# Braker

# Thu Nov 22 09:40:32 2018: Found environment variable $AUGUSTUS_CONFIG_PATH. Setting $AUGUSTUS_CONFIG_PATH to /shelf/augustus_environment/writeable_augustus_config
# Thu Nov 22 09:40:32 2018: Found environment variable $AUGUSTUS_BIN_PATH. Setting $AUGUSTUS_BIN_PATH to /usr/local/Modules/modulefiles/tools/augustus/3.2.2/bin
#conda activate braker

#module load BRAKER2/2.0.4
# THESE TWO HAVE TO GO TOGETHER --hints=hints.gff --optCfgFile=/home/peter/Desktop/R.padi_Braker/extrinsic.bug.cfg

#perl /shelf/apps/User_name/apps/gm_et_linux_64/gmes_petap/gmes_petap.pl --verbose --sequence=/shelf/apps/User_name/newton/RNAseq/braker/GPAL_SCAF_fil_v1UTR_test03/genome.fa --ET=/shelf/apps/User_name/newton/RNAseq/braker/GPAL_SCAF_fil_v1UTR_test03/genemark_hintsfile.gff --cores=16 --soft_mask 1000 1>/shelf/apps/User_name/newton/RNAseq/braker/GPAL_SCAF_fil_v1UTR_test03/GeneMark-ET.stdout 2>/shelf/apps/User_name/newton/RNAseq/braker/GPAL_SCAF_fil_v1UTR_test03/errors/GeneMark-ET.stderr

# without UTR
#/shelf/apps/User_name/conda/envs/braker/bin/perl /shelf/apps/User_name/apps/BRAKER_v2.1.0/braker.pl \

echo "start"
cmd='perl /shelf/apps/User_name/apps/BRAKER_v2.1.0/braker.pl 
--genome=/shelf/apps/User_name/newton/RNAseq/genome.fasta 
--overwrite 
--workingdir=/shelf/apps/User_name/newton/RNAseq/baker_no_UTR
--hints=/shelf/apps/User_name/newton/RNAseq/hintsfile.gff 
--AUGUSTUS_BIN_PATH=/storage/home/users/User_name/shelf_apps/apps/Augustus-master/bin 
--AUGUSTUS_CONFIG_PATH=/storage/home/users/User_name/shelf_apps/apps/Augustus-master/config 
--AUGUSTUS_SCRIPTS_PATH=/storage/home/users/User_name/shelf_apps/apps/Augustus-master/scripts 
--species=GPAL_SCAF_fil_v1_test_004hintsmarvin1 
--augustus_args="--protein=on --start=on --stop=on --cds=on --introns=on --noInFrameStop=true --genemodel=complete " 
--crf --cores 4 --gff3 --rounds 1 --UTR=off --filterOutShort  
--GENEMARK_PATH=/shelf/apps/User_name/apps/gm_et_linux_64/gmes_petap/ '
#echo ${cmd}
#eval ${cmd}


export PATH=/shelf/apps/User_name/apps/BRAKER-master/scripts/:$PATH
# with UTR - doesnt work
echo "start UTRRRRRRRRRR     tests     "
cmd='perl /shelf/apps/User_name/apps/BRAKER-master/scripts/braker.pl
--genome=/shelf/apps/User_name/newton/RNAseq/genome.fasta 
--overwrite 
--workingdir=/shelf/apps/User_name/newton/RNAseq/braker_UTR003
--bam=/shelf/apps/User_name/newton/RNAseq/GpbrakerAligned.bam 
--AUGUSTUS_BIN_PATH=/storage/home/users/User_name/shelf_apps/apps/Augustus-master/bin 
--AUGUSTUS_CONFIG_PATH=/storage/home/users/User_name/shelf_apps/apps/Augustus-master/config 
--AUGUSTUS_SCRIPTS_PATH=/storage/home/users/User_name/shelf_apps/apps/Augustus-master/scripts 
--species=GPAL_SCAF_fil_v1UTR_test0006marvin005
--augustus_args="--protein=on --start=on --stop=on --cds=on --introns=on --noInFrameStop=true --genemodel=complete " 
--crf --cores 4 --gff3 --rounds 1 --UTR=on --filterOutShort  
--GENEMARK_PATH=/shelf/apps/User_name/apps/gm_et_linux_64/gmes_petap/ 
--softmasking '
echo ${cmd}
eval ${cmd}


# with UTR - doesnt work
echo "start gm test"
cmd='perl /shelf/apps/User_name/apps/BRAKER_v2.1.0/braker_PT20181203.pl
--genome=/shelf/apps/User_name/newton/RNAseq/genome.fasta 
--overwrite 
--workingdir=/shelf/apps/User_name/newton/RNAseq/braker_UTR003_Genemarktests
--bam=/shelf/apps/User_name/newton/RNAseq/GpbrakerAligned.bam 
--AUGUSTUS_BIN_PATH=/storage/home/users/User_name/shelf_apps/apps/Augustus-master/bin 
--AUGUSTUS_CONFIG_PATH=/storage/home/users/User_name/shelf_apps/apps/Augustus-master/config 
--AUGUSTUS_SCRIPTS_PATH=/storage/home/users/User_name/shelf_apps/apps/Augustus-master/scripts 
--species=GPAL_SCAF_fil_v1UTR_test0006marvin003gmtest2
--augustus_args="--protein=on --start=on --stop=on --cds=on --introns=on --noInFrameStop=true --genemodel=complete " 
--crf --cores 4 --gff3 --rounds 1 --UTR=on --filterOutShort  
--GENEMARK_PATH=/shelf/apps/User_name/apps/gm_et_linux_64/gmes_petap/ 
--softmasking '
#echo ${cmd}
#eval ${cmd}
