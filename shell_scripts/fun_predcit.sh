# funannotate method
# data must be in ~/programs/funannotate . the following command binds this to the "cd" below.
# https://funannotate.readthedocs.io/en/latest/commands.html
#docker run -it --rm -v $PWD:/home/linuxbrew/data funannotate



####################################################################

cd /home/linuxbrew/data

###########################
# variables
genome="Gp_Newton_haplotype1_sotmasked.fasta"
threads=8
left="R1.fq.gz"
right="R2.fq.gz"
species="Gpal_newtontest03" # two words
strain="newton" # ,ust put a strain
out_folder="Gp_preditc_trinity_braker_pre_run"
max_intron_len=10000
RM_lib="consensi.fa.classified"


################################
#GPAL_newton_softmasked

# 1) # rename the scafold and sort by length
#funannotate clean -i ${genome} --minlen 1000 -o ${genome}.cleaned.fa  

# 2) sort by length
#funannotate sort -i${genome}.cleaned.fa -b scaffold -o ${genome}.cleaned.sorted.fa  

# 3) mask
#funannotate mask -i ${genome}.cleaned.sorted.fa --cpus ${threads} -o MyAssembly.fa  --repeatmodeler_lib ${RM_lib}

 
# 4) train. run trinity and pasa. # if you have stranded --stranded RF,  tight gene density --jaccard_clip,  --nanopore_mrna nanopore_mRNA.fq.gz
# no enough ram on this machine to run trinity
train="funannotate train -i ${genome} -o ${out_folder} --left ${left}  --right ${right} \
	--trinity Trinity-GG.fasta --species ${species} --cpus ${threads} \
	--strain ${strain} --max_intronlen ${max_intron_len}"
#echo ${train}
#eval ${train}



# 5) predict genes - it will know where to look for the out files form the other runs. --organism other (for non fungal)
predict="funannotate predict -i ${genome} -o ${out_folder} \
    --species  ${species} --strain ${strain} --cpus ${threads} --max_intronlen ${max_intron_len} --organism other \
	--busco_db nematoda --database /home/linuxbrew/data/databases/\
	--organism other --augustus_species GPAL_newton_softmasked \
	--transcript_alignments Gpal_newtontest03_newton.valid_blat_alignments.gff3 \
	--rna_bam Gp_soft_maskedAligned.sortedByCoord.out.bam "
echo ${predict} 
eval ${predict}
	
	# possible run with this option --repeats2evm

# 6) if have RNAseq update with get the UTR predctions.
funannotate update -i ${out_folder} --cpus ${threads}


-------------------------------------------------------
Total Gene Models:	19,010
Total transcripts:	19,835
New Gene Models:	99
No Change:		12,243
Update UTRs:		6,604
Exons Changed:		64
Exons/CDS Changed:	0
Dropped Models:		0
CDS AED:		0.006
mRNA AED:		0.043
-------------------------------------------------------
[03:02 PM]: Funannotate update is finished, output files are in the fun_newton/update_results folder
[03:02 PM]: Your next step might be functional annotation, suggested commands:
-------------------------------------------------------
Run InterProScan (Docker required): 
funannotate iprscan -i fun_newton -m docker -c 6

Run antiSMASH: 
funannotate remote -i fun_newton -m antismash -e youremail@server.edu

Annotate Genome: 
funannotate annotate -i fun_newton --cpus 6 --sbt yourSBTfile.txt
-------------------------------------------------------
 

outside of docker:
../funannotate-master/funannotate remote -i fun_newton/ -m all -e User_name@st-andrews.ac.uk --force

run phobius o the server paste result in

back in docker image
           
 funannotate annotate -i fun_newton --cpus 6  --eggnog /home/linuxbrew/data/fun_newton/annotate_misc/genome.proteins.fasta.emapper.annotations --phobius /home/linuxbrew/data/fun_newton/annotate_misc/phobius.results.txt --iprscan /home/linuxbrew/data/fun_newton/annotate_misc/iprscan.xml --busco_db nematoda

####
#  now to annotate
# 7) Interproscan5
#funannotate iprscan -i ${out_folder} -m docker --cpus ${threads}
# may need to run remotely
funannotate remote -i ${out_folder}  -m interproscan -e User_name@st-andrews.ac.uk

# 8) eggnog. Run over net
#http://eggnogdb.embl.de/#/app/emapper 

# 8) final annotate
funannotate annotate -i  ${out_folder} --cpus ${threads}



##################################################################################################
################################

cd /home/linuxbrew/data/funannotate

###########################
# variables
genome="Gp_Newton_haplotype1_sotmasked.fasta"
threads=6
left="R1.fq.gz"
right="R2.fq.gz"
species="Gpal_newtontest03" # two words
strain="newton" # ,ust put a strain
out_folder="fun_${strain}"
max_intron_len=10000
RM_lib="consensi.fa.classified"


################################
#GPAL_newton_softmasked

# 1) # rename the scafold and sort by length
#funannotate clean -i ${genome} --minlen 1000 -o ${genome}.cleaned.fa  

# 2) sort by length
#funannotate sort -i${genome}.cleaned.fa -b scaffold -o ${genome}.cleaned.sorted.fa  

# 3) mask
#funannotate mask -i ${genome}.cleaned.sorted.fa --cpus ${threads} -o MyAssembly.fa  --repeatmodeler_lib ${RM_lib}

 
# 4) train. run trinity and pasa. # if you have stranded --stranded RF,  tight gene density --jaccard_clip,  --nanopore_mrna nanopore_mRNA.fq.gz
# no enough ram on this machine to run trinity
train="funannotate train -i ${genome} -o ${out_folder} --left ${left}  --right ${right} \
	--trinity Trinity-GG.fasta --species ${species} --cpus ${threads} \
	--strain ${strain} --max_intronlen ${max_intron_len}"
echo ${train}
eval ${train}


# 5) predict genes - it will know where to look for the out files form the other runs. --organism other (for non fungal)
predict="funannotate predict -i ${genome} -o ${out_folder} 
    --species  ${species} --strain ${strain} --cpus ${threads} --max_intronlen ${max_intron_len} --organism other 
	--optimize_augustus --augustus_species GPAL_FUNV04f --name GPALNEWT
	--busco_db nematoda --database /home/linuxbrew/data/funannotate/databases/ 
	--organism other --rna_bam Gp_soft_maskedAligned.sortedByCoord.out.bam --other_gff augustus.hints.gff3
	--protein_alignments swiss_prot_align.gff "
echo ${predict}
eval ${predict}

funannotate update -i fun_newton --cpus 6



	# possible run with this option --repeats2evm

# 6) if have RNAseq update with get the UTR predctions.
funannotate update -i ${out_folder} --cpus ${threads}

####
#  now to annotate
# 7) Interproscan5
#funannotate iprscan -i ${out_folder} -m docker --cpus ${threads}
# may need to run remotely
funannotate remote -i ${out_folder}  -m interproscan -e User_name@st-andrews.ac.uk

# 8) eggnog. Run over net
#http://eggnogdb.embl.de/#/app/emapper 

# 8) final annotate
funannotate annotate -i  ${out_folder} --cpus ${threads}






	






	
