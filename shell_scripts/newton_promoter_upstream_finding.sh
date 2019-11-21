
cd /storage/home/users/pjt6/newton/final_genes/upstream/



values=" 300 500 1000 1500 2000 2500 3000 3500 4000 5000 7000 8000 9000 10000 "

cmd="python 
    ~/intergenic_regions/intergenic_regions.py
    --gff ../Gpal_newton_newton.gff3
    -g ../Gpal_newton_newton.scaffolds.fa
    -u 20 
    -z 6
    -o test_20_upstream.fasta "
echo ${cmd}
eval ${cmd}


for v in ${values}:
do
    cmd="python 
    ~/intergenic_regions/intergenic_regions.py
    --gff ../Gpal_newton_newton.gff3
    -g ../Gpal_newton_newton.scaffolds.fa
    -u ${v} 
    -z 200
    -o newton_${v}_upstream.fasta "
    echo ${cmd}
    #eval ${cmd}
done

for v in ${values}:
do
    cmd="python 
    ~/intergenic_regions/intergenic_regions.py
    --gff ../Gpal_newton_newton.gff3
    -g ../Gpal_newton_newton.scaffolds.fa
    -u ${v} 
    -o newton_${v}_upstream.fasta "
    echo ${cmd}
    eval ${cmd}
done

for v in ${values}:
do
    cmd="python 
    ~/intergenic_regions/intergenic_regions.py
    --gff ../Gpal_newton_newton.gff3
    -g ../Gpal_newton_newton.scaffolds.fa
    -u ${v}
    -z 3    
    -o newton_${v}_upstream.fasta "
    echo ${cmd}
    #eval ${cmd}
done


for v in ${values}:
do
    cmd="python 
    get_RBBH_GPPAL_seq.py 
    newton_${v}_upstream.fasta subventral.txt Gp_SUBVENTRAL_${v}_upstream.fasta"
    echo ${cmd}
    eval ${cmd}
done


for v in ${values}:
do
    cmd="python 
    get_sequences_i_want_from_fasta_command_line.py 
    newton_${v}_upstream.fasta GP_non_effectors.names GP_non_effectors_${v}_upstream.fasta"
    echo ${cmd}
    eval ${cmd}
done


# motif finding. 
#FASTA example: findMotifs.pl   targets.fa  fasta   motifResults/   -fasta   background.fa

for v in ${values}:
do
    cmd="findMotifs.pl Gp_SUBVENTRAL_${v}_upstream.fasta fasta homer_results_${v}/ -fasta GP_non_effectors_${v}_upstream.fasta"
    echo ${cmd}
    eval ${cmd}
done

# scan motifs 

#scanMotifGenomeWide.pl motif2.motif  Gp.all_5000_upstream.fasta > motif1_500bp_upstream.bed

