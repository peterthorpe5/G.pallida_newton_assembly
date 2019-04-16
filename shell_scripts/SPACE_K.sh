#!/bin/bash
#$ -cwd
cd /shelf/apps/User_name/newton/newton_wtdgb/L4000

cd /storage/home/users/User_name/shelf_apps/newton/newton_wtdgb/L4000
# purger, purged, output all data, blobtools

#python ~/misc_python/convert_file_format/convert_fq_to_fa.py -i raw_reads.fastq -o raw_reads.fasta

##conda activate python27

#python /shelf/apps/User_name/apps/finishingTool/finisherSC.py -par 16 /storage/home/users/User_name/shelf_apps/newton/newton_wtdgb/L4000/ /shelf/apps/User_name/conda/envs/python27/bin/





#pigz -d ../newton.fastq.gz
perl /shelf/apps/User_name/apps/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE-LongRead.pl -c curated_A80.fasta -p /shelf/apps/User_name/newton/Gp_newton_reduce_haplotypes_400X.correctedReads.fasta -b sspace_long_read_A80_k -t 16 -k 1

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
cd /shelf/apps/User_name/newton/newton_wtdgb/L4000/sspace_long_read_A80_k

/shelf/apps/User_name/apps/gapFinisher/gapFinisher -i /shelf/apps/User_name/newton/newton_wtdgb/L4000/sspace_long_read_A80_k -l /shelf/apps/User_name/newton/Gp_newton_reduce_haplotypes_400X.correctedReads.fasta -m /shelf/apps/User_name/apps/mcr/v717/ -t 1
