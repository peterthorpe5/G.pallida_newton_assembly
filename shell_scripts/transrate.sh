#!/bin/bash
#$ -cwd
#Abort on any error,
#set -e

#echo Running on $HOSTNAME
#echo Current PATH is $PATH
#source $HOME/.bash_profile

################################################################
# Variables: FILLL IN DOWN TO THE END OF VARIABLES

nem_dir=/storage/home/users/User_name/newton/gene_models/fun/final_braker_utr_input_funanot

known_fa="/storage/home/users/User_name/newton/gene_models/fun/Gpal_aa_fromNCBI_JJ.fasta"
known_fa_nucl="Gpal_newton_newton.transcripts.fa"
prefix="GPALNEWT"
test_fa="Gpal_newton_newton.proteins.fa"
test_fa_nt="Gpal_newton_newton.transcripts.fa"

min_len_gene="5"
threads=20
python_directory=$HOME/public_scripts/gene_model_testing/
Working_directory=/storage/home/users/User_name/newton/gene_models/fun/final_braker_utr_input_funanot/
test_gff="${nem_dir}/Gpal_newton_newton.gff3"

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

left_reads="/storage/home/users/User_name/newton/RNAseq/R1.fq.gz"
right_reads="/storage/home/users/User_name/newton/RNAseq/R2.fq.gz"

# END OF VARIABLES !!!!!!!!!!!!!!!!!!!!!!
########################################################
#conda activate genemodel

cd ${Working_directory}



#################################################################################
# change back to the working directory
# Run transrate?
cd ${Working_directory}
echo "Running transrate to get RNAseq mapping and it will tell you 
statistically what genes may be fusions. "
tran="transrate 
	 --assembly Gpal_newton_newton.transcripts.fa 
	 --left ${left_reads} 
	 --right ${right_reads}"
echo ${tran}
eval ${tran}

echo ................................................................ im done
