#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
#$ -N newtool ## job name
#$ -S /bin/bash ## shell where it will run this job
####$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory
###$ -q marvin.q ## queue name
#$ -pe multi 8
#$ -m e
#$ -M User_name@st-andrews.ac.uk

cd /storage/home/users/User_name/shelf_apps/newton

module load oraclejava/jdk1.8.0_74
#default

#/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -i Gp_newton_reduce_haplotypes_400X.correctedReads.fasta.gz -fo prefix

#/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i prefix.ctg.lay -fo prefix.ctg.lay.fa

#/shelf/apps/User_name/canu-1.7.1/Linux-amd64/bin/canu correctedErrorRate=0.85 corMhapSensitivity=normal -p newton_using_EC_reads_EC85 -d newton_using_EC_reads_EC85 genomeSize=120m -pacbio-raw 'Gp_newton_reduce_haplotypes_400X.correctedReads.fasta.gz' -useGrid=False -maxMemory=450 -maxThreads=16


#correctedErrorRate=0.085 corMhapSensitivity=normal


#mv ~/*fastq.gz ./


#/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -i all_pacbio.fastq.gz -fo Hs_raw

#/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i Hs_raw.ctg.lay -fo Hs_raw.ctg.lay.fa


#/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -i newton_all_PacBio.fastq.gz -fo newton_raw

#/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i newton_raw.ctg.lay -fo newton_raw.ctg.lay.fa

#/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -i Hsac_all_data_reduce_ploidy_300X_CMAPSEN_normal.trimmedReads.fasta.gz -fo Hs_trimmed

#/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i Hs_trimmed.ctg.lay -fo Hs_trimmed.ctg.lay.fa


/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -i newton_using_EC_reads_EC85.correctedReads.fasta.gz -fo newt_correc_corre

/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i newt_correc_corre.ctg.lay -fo newt_correc_corre.ctg.lay.fa




