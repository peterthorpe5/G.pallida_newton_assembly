#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
###$ -N canu_16 ## job name
#$ -S /bin/bash ## shell where it will run this job
#$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory
#####$ -q highmemory.q ## queue name
#$ -pe multi 16
#$ -m e
#$ -M User_name@st-andrews.ac.uk
#$ -t 132-140
#$ -tc 11

#cd $HOME/newton/newton_using_EC_reads/correction/1-overlapper
cd /shelf/apps/User_name/newton/newton_using_EC_reads/correction/1-overlapper
#cd /storage/home/users/User_name/newton/newton_using_EC_reads/correction/1-overlapper
#i=$SGE_TASK_ID
#echo $SGE_TASK_ID
./hs_correct_$SGE_TASK_ID.sh
