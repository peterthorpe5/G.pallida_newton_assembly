#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
###$ -N canu_16 ## job name
#$ -S /bin/bash ## shell where it will run this job
#$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory
#####$ -q highmemory.q ## queue name
#$ -pe multi 16
#$ -m e
#$ -M $USER@st-andrews.ac.uk
#$ -t 81-85 
#$ -tc 8

#cd $HOME/newton/newton_using_EC_reads/correction/1-overlapper

cd /storage/home/users/$USER/newton/newton_using_EC_reads/correction/1-overlapper
i=$SGE_TASK_ID
echo $SGE_TASK_ID
./hs_correct_$SGE_TASK_ID.sh
