cd /storage/home/users/pjt6/covid19/Vero-Infected

#java -jar /shelf/training/Trimmomatic-0.38/trimmomatic-0.38.jar PE -summary trim_summary.txt  -threads 16 -phred33 r1.fq.gz r2.fq.gz R1.fq.gz crap1.fastq.gz R2.fq.gz crap2.fastq.gz ILLUMINACLIP:/shelf/training/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:59

perl ~/shelf_apps/apps/trinityrnaseq-Trinity-v2.8.4/util/insilico_read_normalization.pl --left R1.fq --right r2.fq--max_cov 100 --JM 210G --seqType fq--PARALLEL_STATS --cpu 32 --pairs_together
