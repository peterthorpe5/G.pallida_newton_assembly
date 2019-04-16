#!/bin/bash
#$ -cwd

cd /storage/home/users/User_name/shelf_apps/newton/newton_wtdgb/L4000

###################################################################################
#round1

#perl /shelf/apps/User_name/apps/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE-LongRead.pl -g 500 -l 10 -o 100 -c curated_A80.fasta -p /shelf/apps/User_name/newton/Gp_newton_reduce_haplotypes_400X.correctedReads.fasta -b sspace_long_read_A80_k_minlink10_g500_o100 -t 16 -k 1

#cd /shelf/apps/User_name/newton/newton_wtdgb/L4000/sspace_long_read_A80_k_minlink10_g500_o100
#/shelf/apps/User_name/apps/gapFinisher/gapFinisher -i /shelf/apps/User_name/newton/newton_wtdgb/L4000/sspace_long_read_A80_k_minlink10_g500_o100 -l /shelf/apps/User_name/newton/Gp_newton_reduce_haplotypes_400X.correctedReads.fasta -m /shelf/apps/User_name/apps/mcr/v717/ -t 1


###################################################################################
#round2
#perl /shelf/apps/User_name/apps/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE-LongRead.pl -g 500 -l 10 -o 100 -c scaffolds_gapfilled_FINAL.fasta -p /shelf/apps/User_name/newton/Gp_newton_reduce_haplotypes_400X.correctedReads.fasta -b sspace_long_read_A80_k_minlink10_g500_o100 -t 16 -k 1

#cd /shelf/apps/User_name/newton/newton_wtdgb/L4000/sspace_long_read_A80_k_minlink10_g500_o100/sspace_long_read_A80_k_minlink10_g500_o100
#/shelf/apps/User_name/apps/gapFinisher/gapFinisher -i /shelf/apps/User_name/newton/newton_wtdgb/L4000/sspace_long_read_A80_k_minlink10_g500_o100/sspace_long_read_A80_k_minlink10_g500_o100 -l /shelf/apps/User_name/newton/Gp_newton_reduce_haplotypes_400X.correctedReads.fasta -m /shelf/apps/User_name/apps/mcr/v717/ -t 1

#cd /shelf/apps/User_name/newton/newton_wtdgb/L4000/sspace_long_read_A80_k_minlink10_g500_o100/sspace_long_read_A80_k_minlink10_g500_o100

###################################################################################
#round3

#perl /shelf/apps/User_name/apps/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE-LongRead.pl -g 500 -l 10 -o 100 -c scaffolds_gapfilled_FINAL.fasta -p \
#/shelf/apps/User_name/newton/raw_reads.fasta -b sspace_long_read_A80_k_minlink10_g500_o100 -t 16 -k 1

#cd ./sspace_long_read_A80_k_minlink10_g500_o100
#/shelf/apps/User_name/apps/gapFinisher/gapFinisher -i /shelf/apps/User_name/newton/newton_wtdgb/L4000/sspace_long_read_A80_k_minlink10_g500_o100/sspace_long_read_A80_k_minlink10_g500_o100/sspace_long_read_A80_k_minlink10_g500_o100 -l /shelf/apps/User_name/newton/raw_reads.fasta -m /shelf/apps/User_name/apps/mcr/v717/ -t 1

wait
#conda activate purge_haplotigs

cd /shelf/apps/User_name/newton/newton_wtdgb/L4000/sspace_long_read_A80_k_minlink10_g500_o100/sspace_long_read_A80_k_minlink10_g500_o100/sspace_long_read_A80_k_minlink10_g500_o100

minialign -t16 -xont scaffolds_gapfilled_FINAL.fasta /storage/home/users/User_name/shelf_apps/newton/reads/newton_using_EC_reads_EC85.correctedReads.fasta.gz > temp.sam
wait 
samtools view -@ 16 -S -b -o unsorted.bam temp.sam
wait 
samtools sort -@ 16 -o sorted.bam unsorted.bam 

samtools index sorted*bam

#rm temp.sam unsorted.bam

wait 
purge_haplotigs readhist -b sorted*bam -t 8 -g scaffolds_gapfilled_FINAL.fasta
#
#
purge_haplotigs  contigcov  -i sorted*bam.gencov -o coverage_stats.csv  -l 3  -m 27  -h 130
#
##
#
purge_haplotigs purge  -g scaffolds_gapfilled_FINAL.fasta  -c coverage_stats.csv  -b sorted*bam  -t 16 -a 80

perl ~/scaffold_stats.pl -f curated.fasta > curated.stat
