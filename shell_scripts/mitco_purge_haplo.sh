#!/bin/bash
#$ -cwd
cd /shelf/apps/User_name/newton/old_versions/mitoch

conda activate purge_haplotigs


minialign -t16 -xont Gp_Newton_mitochindria.fasta /storage/home/users/User_name/shelf_apps/newton/reads/newton_using_EC_reads_EC85.correctedReads.fasta.gz > temp.sam
wait 
samtools view -@ 16 -S -b -o unsorted.bam temp.sam
wait 
samtools sort -@ 16 -o sorted unsorted.bam 

samtools index sorted*

wait 
purge_haplotigs readhist -b sorted.bam -t 16 -g Gp_Newton_mitochindria.fasta
#
#
purge_haplotigs  contigcov  -i sorted.bam.gencov -o coverage_stats.csv  -l 3  -m 27  -h 130
#
##
#
purge_haplotigs purge  -g Gp_Newton_mitochindria.fasta  -c coverage_stats.csv  -b sorted.bam  -t 16 -a 80

perl ~/scaffold_stats.pl -f curated.fasta | cut -f2

#
#cp curated.fasta curated_A80.fasta
#
#purge_haplotigs purge  -g Gp_Newton_mitochindria.fasta  -c coverage_stats.csv  -b sorted.bam  -t 16 -a 70
#
#cp curated.fasta curated_A70.fasta



#pigz -d ../newton.fastq.gz
#perl /shelf/apps/User_name/apps/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE-LongRead.pl -c curated_A70.fasta -p /shelf/apps/User_name/newton/newton_all_PacBio.fastq -b sspace_long_read_k_A70 -t 16 -k 1

#samtools merge newton.bam \
#/mnt/shared/projects/nematodes/Globodera_pallida/20180507_Gp_Newton/raw_data/r54047_20180424_024920/2_B01/m54047_180424_232338.subreads.bam \
#/mnt/shared/projects/nematodes/Globodera_pallida/20180507_Gp_Newton/raw_data/r54047_20180426_085636/2_B01/m54047_180427_053850.subreads.bam \
#/mnt/shared/projects/nematodes/Globodera_pallida/20180507_Gp_Newton/raw_data/r54047_20180426_085636/3_C01/m54047_180428_015039.subreads.bam \
#/mnt/shared/projects/nematodes/Globodera_pallida/20180507_Gp_Newton/raw_data/r54047_20180426_085636/4_D01/m54047_180428_220428.subreads.bam


#pbbamify --input=sorted_mini2.bam --output=sorted_mini2.pb.bam curated.fasta rawreads.subreads.bam
#cd /home/pt40963/scratch/nematode/newton/
#samtools sort -@ 4 newton.subreads.all.bam sub.sorted
#
#samtools index sub.sorted.bam
#cd /home/pt40963/scratch/nematode/newton/finisher_ALL_haplo_purgedX2_Blobs
#
#pbbamify --input=sorted_mini2.bam --output=sorted_mini.pb.bam curated.fasta sub.sorted.bam
#
#cd /shelf/apps/User_name/newton/newton_wtdgb/L4000/sspace_long_read_k

#/shelf/apps/User_name/apps/gapFinisher/gapFinisher -i /shelf/apps/User_name/newton/newton_wtdgb/L4000/sspace_long_read_k_A70 -l /shelf/apps/User_name/newton/newton_wtdgb/L4000/raw_reads.fasta  -m /shelf/apps/User_name/apps/mcr/v717/ -t 1
