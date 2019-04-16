#!/bin/bash
#$ -cwd
#set -e

cd /shelf/apps/User_name/newton/stricter_scaff_polish

conda activate python27

# HTseq counts
cd /storage/home/users/User_name/newton/final_genome2/RNAseq_mapping/

#conda activate python27
#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_14DPI_rep1*.bam Gpal_newton_newton.gff3 > Gp_14DPI_rep1_genes.counts
#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_14DPI_rep2*.bam Gpal_newton_newton.gff3 > Gp_14DPI_rep2_genes.counts
#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_21DPI_rep1*.bam Gpal_newton_newton.gff3 > Gp_21DPI_rep1_genes.counts


#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_21DPI_rep2*.bam Gpal_newton_newton.gff3 > Gp_21DPI_rep2_genes.counts
htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_28DPI_rep1*.bam Gpal_newton_newton.gff3 > Gp_28DPI_rep1_genes.counts

#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_28DPI_rep2*.bam Gpal_newton_newton.gff3 > Gp_28DPI_rep2_genes.counts
#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_35DPI_rep1*.bam Gpal_newton_newton.gff3 > Gp_35DPI_rep1_genes.counts

#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_35DPI_rep2*.bam Gpal_newton_newton.gff3 > Gp_35DPI_rep2_genes.counts
#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_7DPI_rep1*.bam Gpal_newton_newton.gff3 > Gp_7DPI_rep1_genes.counts

#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_7DPI_rep2*.bam Gpal_newton_newton.gff3 > Gp_7DPI_rep2_genes.counts
#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_EGG_rep1*.bam Gpal_newton_newton.gff3 > Gp_EGG_rep1_2_genes.counts

#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_EGG_rep2*.bam Gpal_newton_newton.gff3 > Gp_EGG_rep2_genes.counts
#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_J2_rep1*.bam Gpal_newton_newton.gff3 > Gp_J2_rep1_genes.counts
#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_J2_rep2*.bam Gpal_newton_newton.gff3 > Gp_J2_rep2_genes.counts

#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_J2_rep3*.bam Gpal_newton_newton.gff3 > Gp_J2_rep3_genes.counts
#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_MALE_rep1*.bam Gpal_newton_newton.gff3 > Gp_MALE_rep1_genes.counts
#htseq-count  -s no  -r pos --quiet -f bam --type gene --idattr ID Gp_MALE_rep2*.bam Gpal_newton_newton.gff3 > Gp_MALE_rep2_1_genes.counts

echo "	Gp_7DPI_rep1	Gp_7DPI_rep2	Gp_14DPI_rep1	Gp_14DPI_rep2	Gp_21DPI_rep1	Gp_21DPI_rep2	Gp_28DPI_rep1	Gp_28DPI_rep2	Gp_35DPI_rep1	Gp_35DPI_rep2	Gp_EGG_rep1	Gp_EGG_rep2	Gp_J2_rep1	Gp_J2_rep2	Gp_J2_rep3	Gp_MALE_rep1	Gp_MALE_rep2" > Gp_genes.counts.matrix
# merge these into a file
FILES=$(ls -t -v *.counts | tr '\n' ' ')
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES | grep -v "__" >> Gp_genes.counts.matrix
