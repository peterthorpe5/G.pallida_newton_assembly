#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
#$ -N canu_16 ## job name
#$ -S /bin/bash ## shell where it will run this job
####$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory
#$ -q highmemory.q ## queue name
#$ -pe multi 16
#$ -m e
#$ -M pjt6@st-andrews.ac.uk

cd $HOME/newton

module load oraclejava/jdk1.8.0_74
#default

/storage/home/users/pjt6/canu/Linux-amd64/bin/canu -p newton_using_EC_reads correctedErrorRate=0.15 cnsErrorRate=0.25 -d newton_using_EC_reads genomeSize=120m -pacbio-raw 'Gp_newton_reduce_haplotypes_400X.correctedReads.fasta.gz' -useGrid=False -maxMemory=450 -maxThreads=16
