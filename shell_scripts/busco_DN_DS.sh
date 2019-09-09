#!/usr/bin/env bash
#! -cwd
#
# run_codon_phyml_dn_ds.sh
#
# shell to run codonphyml on all the files in a folder.
#

#set -e
# Define input files and input/output directories

cd /storage/home/users/pjt6/newton/comparative_genomics/all_nems/busco/
rm *fasta_*
rm -rf back_translated
filenames=*.fasta



for f in ${filenames}
do
	#echo "Running muscle ${f}"
	cmd="muscle -in ${f} -out ${f}_aligned.fasta" 
	echo ${cmd}
	eval ${cmd}
	wait
	
done
###############################
###############################
# purge bad aligning sequences: https://github.com/dukecomeback/bad-sequence-remover


filenames2=*aligned.fasta
for file in ${filenames2}
do
	echo "Running seq purging ${f}"
	cmd="perl /storage/home/users/pjt6/bad-sequence-remover/badseq_remover.pl
    -in ${file} 
    -out ${file}_purges.fasta
    -exp GPALN" 
	echo ${cmd}
	eval ${cmd}
	wait
done

###############

filenames2=*_purges.fasta
for file in ${filenames2}
do
	echo "Running muscle ${f}"
	cmd="muscle -in ${file} -out ${file}_refine.fasta -refine" 
	echo ${cmd}
	eval ${cmd}
	wait
done
mkdir back_translated

python /storage/home/users/pjt6/newton/comparative_genomics/all_nems/Align_back_translate.py

#####################
cd back_translated/
#filenames3=*.fasta
#for file in ${filenames3}
#do
#	echo "Running trimal ${file}"
#	cmd="/storage/home/users/pjt6/trimAl/source/trimal  
#    -in ${file}  
#    -out ${file}.phy
#    -phylip 
#    -gappyout" 
#	echo ${cmd}
#	eval ${cmd}
#	wait
#done

cd back_translated/
filenames3=*translated.fasta
for file in ${filenames3}
do
	echo "Running trimal ${file}"
	cmd="/storage/home/users/pjt6/trimAl/source/trimal  
    -in ${file}  
    -out ${file}.phy
    -phylip 
    -nogaps" 
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

grep "Probability: o3=" -H *_stats.txt > ~/newton/comparative_genomics/all_nems/busco_no_gaps_NO_MINC_orthofinder_DN_DS.txt





