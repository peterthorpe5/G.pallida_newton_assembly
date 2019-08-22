#!/usr/bin/env bash
#! -cwd
#
# run_codon_phyml_dn_ds.sh
#
# shell to run codonphyml on all the files in a folder.
#

#set -e
# Define input files and input/output directories

cd /storage/home/users/pjt6/newton/comparative_genomics/all_nems/sub

filenames=*.fasta





for f in ${filenames}
do
	echo "Running muscle ${f}"
	cmd="muscle -in ${f} -out ${f}_aligned.fasta" 
	echo ${cmd}
	eval ${cmd}
	wait
	
done
python ../Align_back_translate_Aug2014_cow_boy_method002.py

filenames2 = *_aligned.fasta
for file in ${filenames2}
do
	echo "Running muscle ${f}"
	cmd="muscle -in ${file} -out ${file}_refine.fasta -refine" 
	echo ${cmd}
	eval ${cmd}
	wait
	
done




