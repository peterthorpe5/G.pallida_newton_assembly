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
threads=60
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


###########################################################################

# --taxonmap /shelf/public/blastntnr/blastDatabases//prot.accession2taxid.gz	
echo "running diamond-BLAST against NR xml for BLAST2GO"
diam_p="diamond blastp 
		-p ${threads} 
		--more-sensitive -e 0.00001 
	   -q aa.fa
	   -d ${diamond_db}  
	   --outfmt 5 --salltitles
	   --out aa2.fasta_vs_nr_B2GO.xml"
echo ${diam_p}
eval ${diam_p}
wait
