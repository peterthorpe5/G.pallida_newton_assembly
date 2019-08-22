#!/usr/bin/env bash
#! -cwd
#
# run_codon_phyml_dn_ds.sh
#
# shell to run codonphyml on all the files in a folder.
#

#set -e
# Define input files and input/output directories

cd /storage/home/users/pjt6/newton/comparative_genomics/all_nems/effectors

filenames=*.fasta



for f in ${filenames}
do
	#echo "Running muscle ${f}"
	cmd="muscle -in ${f} -out ${f}_aligned.fasta" 
	#echo ${cmd}
	#eval ${cmd}
	wait
	
done
###############################
python ../Align_back_translate_Aug2014_cow_boy_method002.py

filenames2 = *_aligned.fasta
for file in ${filenames2}
do
	#echo "Running muscle ${f}"
	cmd="muscle -in ${file} -out ${file}_refine.fasta -refine" 
	#echo ${cmd}
	#eval ${cmd}
	wait
	
done

#####################
cd back_translated/
filenames3=*.fasta
for file in ${filenames3}
do
	echo "Running trimal ${file}"
	cmd="/storage/home/users/pjt6/trimAl/source/trimal  
    -in ${file}  
    -out ${file}.phy
    -phylip 
    -gappyout" 
	echo ${cmd}
	eval ${cmd}
	wait
done

mkdir phy_files

mv *.phy ./phy_files/

cd phy_files


filenames=*.phy

# phyml
export PATH=//storage/home/users/pjt6/phyml-20140723/src/:$PATH

#codonphyml
export PATH=/storage/home/users/pjt6/codonPHYML_dev/src/:$PATH



## nucleotide version
mkdir files_done
for f in ${filenames}
do
	echo "Running codonphyml ${f}"
	cmd="codonphyml -i ${f} -m GY --fmodel F3X4 -t e -f empirical -w g -a e" #--optBrent 3
	echo ${cmd}
	eval ${cmd}
	wait
	cmd2="mv ${f} ./files_done"
	echo ${cmd2}
	eval ${cmd2}
	
done
