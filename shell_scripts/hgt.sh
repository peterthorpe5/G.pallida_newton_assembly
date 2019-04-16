#!/bin/bash
#$ -cwd
#Abort on any error,
#set -e

#echo Running on $HOSTNAME
#echo Current PATH is $PATH
#source $HOME/.bash_profile

################################################################
# Variables: FILLL IN DOWN TO THE END OF VARIABLES

nem_dir=/storage/home/users/User_name/newton/gene_models/fun/incorrect_models_trin_pasa_braker_passed_to_it

known_fa="${nem_dir}/Gpal_aa_fromNCBI_JJ.fasta"
known_fa_nucl="Gpal_newtontest03_newton.transcripts.fa"
prefix="GPALNEWT"
test_fa="Gpal_newtontest03_newton.proteins.fa"
test_fa_nt="Gpal_newtontest03_newton.transcripts.fa"

min_len_gene="5"
threads=8
python_directory=$HOME/public_scripts/gene_model_testing/
Working_directory=/storage/home/users/User_name/newton/gene_models/fun/incorrect_models_trin_pasa_braker_passed_to_it
test_gff="${nem_dir}/Gpal_newtontest03_newton.gff3"

genome="/storage/home/users/User_name/newton/gene_models/Gp_Newton_haplotype1.fasta"

diamond_db=/shelf/public/blastntnr/blastDatabases/nr_PT.dmnd

# FOR HGT
# tax_filter_out is the phylum your beast lives in, or clade if you want to get a more refined HGT result
# tax_filter_up_to e.g. metazoan = tax_filter_up_to
# for aphid: #tax_filter_out=6656 #tax_filter_up_to=33208 
# for nematodes,
tax_filter_out=6231
tax_filter_up_to=33208
species_tx_id=36090 


# If you want to run transrate to get the RNAseq read mapping to gene 
# fill these out. Else, just let it fail, as it is the last step.

left_reads="/storage/home/users/User_name/newton/gene_models/fun/unmapped_R1.fq.gz"
right_reads="/storage/home/users/User_name/newton/gene_models/fun/unmapped_R2.fq.gz"

# END OF VARIABLES !!!!!!!!!!!!!!!!!!!!!!
########################################################
#conda activate genemodel

cd ${Working_directory}


echo "adding tx_id and descriptions to diamond-BLAST output"
tax="python $HOME/public_scripts/Diamond_BLAST_add_taxonomic_info/Diamond_blast_to_taxid.py
	-i aa.fasta_vs_nr.tab 
	-p /shelf/public/blastntnr/blastDatabases/ 
	-o aa.fasta_vs_nr_tax.tab"
#echo ${tax}
#eval ${tax}
wait

echo "predicting HGT"
HGT="python $HOME/public_scripts/Lateral_gene_transfer_prediction_tool/Lateral_gene_transfer_predictor.py 
		-i *_vs_nr_tax.tab 
		--tax_filter_out ${tax_filter_out} 
		--tax_filter_up_to ${tax_filter_up_to}
		-p /shelf/public/blastntnr/blastDatabases/ -o LTG_results.out"

echo ${HGT}
eval ${HGT}
wait
#

#Filter taxomony commands:
echo "filtering blast results"
filter_top_blasts="python ${python_directory}/top_BLAST_hit_filter_out_tax_id.py 
				  -i *_vs_nr_tax.tab 
				  -t ${tax_filter_out} 
				  -p /shelf/public/blastntnr/blastDatabases/ 
				  -o top_not_phylum_${tax_filter_out}.hits"
echo ${filter_top_blasts}
eval ${filter_top_blasts}
wait

filter_species="python ${python_directory}/top_BLAST_hit_filter_out_tax_id.py 
			   -i *_vs_nr_tax.tab 
			   -t ${species_tx_id}
			   -p /shelf/public/blastntnr/blastDatabases/ 
			   -o top_not_species_tx_id_${species_tx_id}.hits"
echo ${filter_species}
eval ${filter_species}
wait


