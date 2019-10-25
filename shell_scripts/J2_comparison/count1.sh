#!/bin/bash
#$ -cwd
#set -e

cd /storage/home/users/pjt6/newton/RNAseq_J2

conda activate python27

# HTseq counts

#conda activate python27
htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID \
Gp*.bam Gpal_newton_newton.gff3 > Gp_J2_newton_genes.counts
