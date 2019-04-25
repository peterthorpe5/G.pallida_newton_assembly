#!/bin/bash
#$ -cwd
#$ -V

conda activate homer

values="500 1000 1500 2000 2500 3000 5000"

for v in ${values}:
do
    motif=1
    #cmd="findMotifs.pl Gp_SUBVENTRAL_${v}_upstream.fasta fasta homer_results__${v}/ -fasta GP_non_effectors_${v}_upstream.fasta"
    awk '{print $1}' Gp.all_${v}_upstream.fasta > Gp.all_${v}_upstream_no_des.fasta
    cmd="scanMotifGenomeWide.pl 
    ./BLAST_hit_candiates_sig_p_upJ2_homer_results_${v}/homerResults/motif${motif}.motif
    Gp.all_${v}_upstream_no_des.fasta 
    > 
    upstream_${v}_motif${motif}_${v}bp_upstream.bed"
    echo ${cmd}
    eval ${cmd}
    pycmd="/shelf/apps/pjt6/conda/envs/python36/bin/python 
    parse_scan_motifs.py 
    -i upstream_${v}_motif${motif}_${v}bp_upstream.bed
    -o upstream_${v}_motif${motif}_${v}bp_upstream.bed.summary"
    echo ${pycmd}
    eval ${pycmd}
done


for v in ${values}:
do
    motif=2
    #cmd="findMotifs.pl Gp_SUBVENTRAL_${v}_upstream.fasta fasta homer_results__${v}/ -fasta GP_non_effectors_${v}_upstream.fasta"
    awk '{print $1}' Gp.all_${v}_upstream.fasta > Gp.all_${v}_upstream_no_des.fasta
    cmd="scanMotifGenomeWide.pl 
    ./BLAST_hit_candiates_sig_p_upJ2_homer_results_${v}/homerResults/motif${motif}.motif
    Gp.all_${v}_upstream_no_des.fasta 
    > 
    upstream_${v}_motif${motif}_${v}bp_upstream.bed"
    echo ${cmd}
    eval ${cmd}
    pycmd="/shelf/apps/pjt6/conda/envs/python36/bin/python 
    parse_scan_motifs.py 
    -i upstream_${v}_motif${motif}_${v}bp_upstream.bed
    -o upstream_${v}_motif${motif}_${v}bp_upstream.bed.summary"
    echo ${pycmd}
    eval ${pycmd}
done

for v in ${values}:
do
    motif=3
    #cmd="findMotifs.pl Gp_SUBVENTRAL_${v}_upstream.fasta fasta homer_results__${v}/ -fasta GP_non_effectors_${v}_upstream.fasta"
    awk '{print $1}' Gp.all_${v}_upstream.fasta > Gp.all_${v}_upstream_no_des.fasta
    cmd="scanMotifGenomeWide.pl 
    ./BLAST_hit_candiates_sig_p_upJ2_homer_results_${v}/homerResults/motif${motif}.motif
    Gp.all_${v}_upstream_no_des.fasta 
    > 
    upstream_${v}_motif${motif}_${v}bp_upstream.bed"
    echo ${cmd}
    eval ${cmd}
    pycmd="/shelf/apps/pjt6/conda/envs/python36/bin/python 
    parse_scan_motifs.py 
    -i upstream_${v}_motif${motif}_${v}bp_upstream.bed
    -o upstream_${v}_motif${motif}_${v}bp_upstream.bed.summary"
    echo ${pycmd}
    eval ${pycmd}
done
