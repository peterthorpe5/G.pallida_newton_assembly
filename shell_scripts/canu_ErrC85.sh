#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
#$ -N canu_EC85 ## job name
#$ -S /bin/bash ## shell where it will run this job
####$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory
#$ -q marvin.q ## queue name
#$ -pe multi 16
#$ -m e
#$ -M User_name@st-andrews.ac.uk

cd /storage/home/users/User_name/shelf_apps/newton

module load oraclejava/jdk1.8.0_74
#default

/shelf/apps/User_name/canu-1.7.1/Linux-amd64/bin/canu correctedErrorRate=0.85 corMhapSensitivity=normal -p newton_using_EC_reads_EC85 -d newton_using_EC_reads_EC85 genomeSize=120m -pacbio-raw 'Gp_newton_reduce_haplotypes_400X.correctedReads.fasta.gz' -useGrid=False -maxMemory=450 -maxThreads=16


#correctedErrorRate=0.085 corMhapSensitivity=normal
