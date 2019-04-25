#!/bin/bash

#$ -cwd
#$ -V
cd /storage/home/users/pjt6/newton/RNAseq
conda activate python36
module load samtools

#cp /storage/home/users/pjt6/newton/final_genes/funaannot/update_results/Gpal_newton_newton.gff3 ./


python /storage/home/users/pjt6/public_scripts/TransStart/TransStart.py -b Gp_finalAligned_not_masked.sortedByCoord.out.bam \
-g Gp_Newton_haplotype1.fasta \
--gff Gpal_newton_newton.gff3 -o Gpal_newton_TranStart.txt \
 --keep_gene_depth yes 


 python /storage/home/users/pjt6/public_scripts/TransStart/GFF_to_fasta.py --gff Gpal_newton_TranStartbased_on_min_value.gff -g Gp_Newton_haplotype1.fasta -o Gpal_newton_transtart.fasta

