#!/bin/bash
#$ -cwd
#Abort on any error,
#set -e

#echo Running on $HOSTNAME
#echo Current PATH is $PATH
#source $HOME/.bash_profile

################################################################
# Variables: FILLL IN DOWN TO THE END OF VARIABLES

nem_dir=/storage/home/users/User_name/newton/gene_models/fun/trinity_pasa_RNAseq/

known_fa="${nem_dir}/Gpal_aa_fromNCBI_JJ.fasta"
known_fa_nucl="Gpal_newtontest03_newton.transcripts.fa"
prefix="GPALNEWT"
test_fa="Gpal_newtontest03_newton.proteins.fa"
test_fa_nt="Gpal_newtontest03_newton.transcripts.fa"

min_len_gene="5"
threads=8
python_directory=$HOME/public_scripts/gene_model_testing/
Working_directory=/storage/home/users/User_name/newton/gene_models/fun/trinity_pasa_RNAseq/
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


# If you want to run transrate to get the RNAseq read mapping to gene 
# fill these out. Else, just let it fail, as it is the last step.

left_reads="/storage/home/users/User_name/newton/gene_models/fun/unmapped_R1.fq.gz"
right_reads="/storage/home/users/User_name/newton/gene_models/fun/unmapped_R2.fq.gz"

# END OF VARIABLES !!!!!!!!!!!!!!!!!!!!!!
########################################################
#conda activate genemodel

cd ${Working_directory}

# for now in testing
 rm -rf gff_stats known_fa_all_hits


echo "For this program you need:
Genome tools. Blast. gffread (cufflinks). Python. Biopython. diamond"

mkdir gff_stats

# genome tools to check and reformat the gff - essential step!
echo "prepare the amino acid seq from the GFF. First check gff"
gt_cmd="gt gff3 -sort -retainids -tidy -addids -sortnum -fixregionboundaries 
	   -addintrons -force -o ${test_gff}_reformatted.gff3 ${test_gff}"
echo ${gt_cmd}
eval ${gt_cmd}
wait

# GT to get stats on the predicted models. 

gt="gt gff3 -sort -tidy -addintrons ${test_gff}_reformatted.gff3 | 
   gt stat -force -genelengthdistri -o augustus_genelengthdistri.STAT 
   > temp"
echo ${gt}
eval ${gt}
wait
gt="gt gff3 -sort -tidy -addintrons ${test_gff}_reformatted.gff3 | 
gt stat -force -genescoredistri -o augustus_genescoredistri.STAT 
> temp"
echo ${gt}
eval ${gt}
wait
gt="gt gff3 -sort -tidy -addintrons ${test_gff}_reformatted.gff3 
| gt stat -force -exonlengthdistri -o augustus_exonlengthdistri.STAT 
> temp"
echo ${gt}
eval ${gt}
wait
gt="gt gff3 -sort -tidy -addintrons ${test_gff}_reformatted.gff3 
| gt stat -force -exonnumberdistri -o augustus_exonnumberdistri.STAT 
> temp"
echo ${gt}
eval ${gt}
wait
gt="gt gff3 -sort -tidy -addintrons ${test_gff}_reformatted.gff3 
| gt stat -force -intronlengthdistri -o augustus_intronlengthdistri.STAT 
> temp"
echo ${gt}
eval ${gt}
wait
gt="gt gff3 -sort -tidy -addintrons ${test_gff}_reformatted.gff3 
| gt stat -force -cdslengthdistri -o augustus_cdslengthdistri.STAT 
> temp"
echo ${gt}
eval ${gt}
wait

mv *.STAT ./gff_stats

rm temp

# rename gff genes
echo "step: rename the gene in the reformatted GFF file"
rename_gff_genes="python ${python_directory}/re_name_genes.py 
		        --gff ${test_gff}_reformatted.gff3 
		        --prefix ${prefix}
		        -o ${prefix}"
echo ${rename_gff_genes}
eval ${rename_gff_genes}
wait

# gffread to the the cds from the reformat the gff - essential step!
echo "getting the AA and nt cds from the genome"
gff_to_cds="gffread ${test_gff} 
		   -g ${genome} -x nt.fa 
		   -y temp.fa"
#echo ${gff_to_cds}
#eval ${gff_to_cds}
wait
cp ${test_fa} temp.fa
cp ${test_fa_nt} nt.fa 
# remove dots as stop codons
echo "step: rename the gene in the reformatted GFF file"
remove_stops="python ${python_directory}/rewrite_as_fasta.py 
		      -i temp.fa
			  -l ${min_len_gene}
		      -o aa.fa"
echo ${remove_stops}
eval ${remove_stops}
wait

# going to change the input gene models to this one without the "*" stops.
#  As this breakseverything. 
test_fa="aa.fa"

# make blastdb
echo "step1: make blastdb"
mkdb="makeblastdb -in ${test_fa} -dbtype prot"
echo ${mkdb}
eval ${mkdb}
wait

#blast to xml
echo "step2: blast to xml"
bl_p="blastp -db ${test_fa} -query ${known_fa} -evalue 1e-10 
	  -seg no -num_threads ${threads} 
      -outfmt 5 -out test_fa_vs_known_fa.xml"
echo ${bl_p}
eval ${bl_p}
wait

mkdir known_fa_all_hits

#blast 2 tab
echo "step2b: blast to tab - top hit only"
bl_p2="blastp -db ${test_fa} -query ${known_fa} -evalue 1e-10 
	  -seg no -max_target_seqs 1 -num_threads ${threads} 
	  -outfmt 7 -out test_fa_vs_known_fa.tab"
echo ${bl_p2}
eval ${bl_p2}
wait

#blast 2 tab
echo "step2b: blast to tab - top hit only"
no_comment="grep -v '#' test_fa_vs_known_fa.tab > test_fa_vs_known_fa2.tab "
echo ${no_comment}
eval ${no_comment}
wait

# graph the blast results
graph="python ${python_directory}/blast_stats.py 
	  -i test_fa_vs_known_fa.tab 
	  -o test_fa_vs_known_fa.graphs"
echo ${graph}
eval ${graph}
wait


# convert the xml
echo "step3: convert the xml file"
parse_xml="python ${python_directory}/BLAST_parser_return_hits_NAME_only.py 
		  -i test_fa_vs_known_fa.xml 
		  -o test_fa_vs_known_fa.xml.condensed.out"
echo ${parse_xml}
eval ${parse_xml}
wait

# get the matches seq for the tab blast to the effectors
echo "step4: get the seqs of the top hit."
get_seq="python ${python_directory}/Get_sequence_from_tab_blast.py 
	    -b test_fa_vs_known_fa.tab 
		--known_fa  Gpal_aa_fromNCBI_JJ.fasta 
		--folder known_fa_all_hits
		-p ${test_fa} 
		-n ${test_fa}"
echo ${get_seq}
eval ${get_seq}
wait

# change to where the file have been put
cd ./known_fa_all_hits

# delete empty files
delete_empty="find . -size 0 -delete"
echo ${delete_empty}
eval ${delete_empty}
wait

filenames=*.fasta

for f in ${filenames}
do
	echo "Running muscle ${f}"
	cmd="muscle3.8.31_i86linux64 -in ${f} 
		-out ${f}_aligned.fasta -maxiters 5000 -maxtrees 15" 
	echo ${cmd}
	eval ${cmd}
	wait
done

filenames2=*_aligned.fasta
for file in ${filenames2}
do
	echo "Running muscle ${f}"
	cmd="muscle3.8.31_i86linux64 
		-in ${file} -out ${file}_refine.fasta 
		-refine" 
	echo ${cmd}
	eval ${cmd}
	wait
done
	
rm *_aligned.fasta

mkdir alignments
mv *_refine.fasta ./alignments

###########################################################################
# diamond blast aginst NR?
cd ${Working_directory}

echo "running diamond-BLAST against NR"
diam_p="diamond blastp -p ${threads} --more-sensitive -e 0.00001 
	   -v -q aa.fa 
	   -d ${diamond_db}  
	   -a aa.fasta_vs_nr.da"
echo ${diam_p}
eval ${diam_p}
wait

# --taxonmap /shelf/public/blastntnr/blastDatabases//prot.accession2taxid.gz	
echo "running diamond-BLAST against NR xml for BLAST2GO"
diam_p="diamond blastp 
		-p ${threads} 
		--more-sensitive -e 0.00001 
	   -q GCA_003473485.2_MtrunA17r5.0-ANR_protein.faa
	   -d ${diamond_db}  
	   --outfmt 5 --salltitles
	   --out aa.fasta_vs_nr.xml"
#echo ${diam_p}
#eval ${diam_p}
wait

echo "converting diamond-BLAST output"
diam_v="diamond view -a aa.fasta*.daa -f tab -o aa.fasta_vs_nr.tab"
echo ${diam_v}
eval ${diam_v}
wait

echo "converting diamond-BLAST output to xml for blast2 go"
diam_v="diamond view -a aa.fasta*.daa --outfmt 5 -o aa.fasta_vs_nr.xml"
echo ${diam_v}
eval ${diam_v}
wait

echo "adding tx_id and descriptions to diamond-BLAST output"
tax="python $HOME/public_scripts/Diamond_BLAST_add_taxonomic_info/Diamond_blast_to_taxid.py
	-i aa.fasta_vs_nr.tab 
	-p /shelf/public/blastntnr/blastDatabases/ 
	-o aa.fasta_vs_nr_tax.tab"
echo ${tax}
eval ${tax}
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


#################################################################################
# change back to the working directory
# Run transrate?
cd ${Working_directory}
echo "Running transrate to get RNAseq mapping and it will tell you 
statistically what genes may be fusions. "
tran="transrate 
	 --assembly nt.fa 
	 --left ${left_reads} 
	 --right ${right_reads}"
echo ${tran}
eval ${tran}

echo ................................................................ im done
