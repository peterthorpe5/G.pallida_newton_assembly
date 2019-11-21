#!/bin/bash
#$ -cwd
#$ -V

conda activate homer

values=" 300 500 1000 1500 2000 2500 3000 5000"
motif=" 1 2 3 "

for v in ${values}:
do
    for j in ${values}:
    do
        for m in ${motif}:
        do
            #cmd="findMotifs.pl Gp_SUBVENTRAL_${v}_upstream.fasta fasta homer_results__${v}/ -fasta GP_non_effectors_${v}_upstream.fasta"
            awk '{print $1}' newton_${j}_upstream.fasta | grep -v "T2" > Gp.all_${j}_upstream_no_des.fasta
            cmd="scanMotifGenomeWide.pl 
            ./homer_results_${v}/homerResults/motif${m}.motif
            Gp.all_${j}_upstream_no_des.fasta
            > 
            upstream_${v}_motif${m}_${j}bp_upstream.bed"
            echo ${cmd}
            eval ${cmd}
            pycmd="/shelf/apps/pjt6/conda/envs/python36/bin/python 
            parse_scan_motifs.py 
            -i upstream_${v}_motif${m}_${j}bp_upstream.bed
            -o upstream_${v}_motif${m}_${j}bp_upstream.bed.summary"
            echo ${pycmd}
            eval ${pycmd}
        done
    done
done

