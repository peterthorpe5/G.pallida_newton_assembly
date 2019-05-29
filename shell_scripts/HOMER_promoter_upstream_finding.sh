#!/bin/bash
#$ -cwd
#$ -V

cd /storage/home/users/pjt6/newton/final_genes/funaannot/update_results/upstream_regions

#conda activate python27
# prepare the gene names file

# prepare the gff file
#python ~/public_scripts/genomic_upstream_regions/re_format_GFF_Mcscanx.py --gff ../Gpal_newton_newton.gff3 -m True -o gp.format.out

values="500 1000 1500 2000 2500 3000 5000 10000"

for v in ${values}:
do
    cmd="python ~/public_scripts/genomic_upstream_regions/get_upstream_regions.py 
    -c gp.format.out 
    -g ../Gpal_newton_newton.scaffolds.fa 
    -f all_gene.names 
    -u ${v} 
    -o Gp.all_${v}_upstream.fasta 
    > warning_${v}.out"
    echo ${cmd}
    eval ${cmd}
done

# motif finding. 


for v in ${values}:
do
    cmd="python 
    get_RBBH_GPPAL_seq.py 
    Gp.all_${v}_upstream.fasta subventral.txt Gp_SUBVENTRAL_${v}_upstream.fasta"
    echo ${cmd}
    eval ${cmd}
done


for v in ${values}:
do
    cmd="python 
    get_sequences_i_want_from_fasta_command_line.py 
    Gp.all_${v}_upstream.fasta GP_non_effectors.names GP_non_effectors_${v}_upstream.fasta"
    echo ${cmd}
    eval ${cmd}
done



conda activate homer

#FASTA example: findMotifs.pl   targets.fa  fasta   motifResults/   -fasta   background.fa

for v in ${values}:
do
    cmd="findMotifs.pl   
    Gp_SUBVENTRAL_${v}_upstream.fasta 
    homer_results_${v}/ 
    -fasta 
    GP_non_effectors_${v}_upstream.fasta"
    echo ${cmd}
    eval ${cmd}
done

#findMotifs.pl   Gp_SUBVENTRAL_1000_upstream.fasta homer_results/ -fasta GP_non_effectors_1000_upstream.fasta

