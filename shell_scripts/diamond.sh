#! /bin/bash $ -cwd Abort on any error, set -e module load diamond
# for nematoda cd /storage/home/users/User_name/medicago/A17
#cd /storage/home/users/User_name/medicago/DZA/gene_predictions
cd /storage/home/users/User_name/newton/gene_models/fun/final_genes/annotations
# --taxonmap $HOME/scratch/blast_databases/prot.accession2taxid.gz
echo "running diamond-BLAST against NR" 
diam_p="/storage/home/users/User_name/shelf_apps/apps/diamond blastp
	-p 8 --more-sensitive -e 0.00001
	 -q GPAL_newton_V1.0.proteins.fa
	   -d /shelf/public/blastntnr/blastDatabases/nr_PT.dmnd
	   --outfmt 5 --salltitles
	   --out GPAL_newton_V1.0.proteins.fa_vs_nr.xml" 
echo ${diam_p} 
eval ${diam_p}
wait
