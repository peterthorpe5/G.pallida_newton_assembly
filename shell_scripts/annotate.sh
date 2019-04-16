# funannotate method
# data must be in ~/programs/funannotate . the following command binds this to the "cd" below.
# https://funannotate.readthedocs.io/en/latest/commands.html
#docker run -it --rm -v $PWD:/home/linuxbrew/data funannotate
set -e


####################################################################
cd /home/User_name/programs/funannotate1.5.2


###########################
# variables
genome="Gp_Newton_haplotype1.softmaked.fasta"
threads=6
left="R1.fq.gz"
right="R2.fq.gz"
species="Gpal_newton" # two words
strain="newton" # ,ust put a strain
out_folder="final_newton2"
max_intron_len=10000
RM_lib="consensi.fa.classified"
bam="Gp_finalAligned.sortedByCoord.out.bam"


 

#outside of docker:

~/programs/fun_master/funannotate/funannotate remote -i final_newton2/ -m interproscan -e ????@st-andrews.ac.uk --force

# 8) eggnog. Run over net
#http://eggnogdb.embl.de/#/app/emapper 
run phobius o the server paste result in

#back in docker image
           
# funannotate annotate -i fun_newton --cpus 6  --eggnog /home/linuxbrew/data/fun_newton/annotate_misc/genome.proteins.fasta.emapper.annotations --phobius /home/linuxbrew/data/fun_newton/annotate_misc/phobius.results.txt --iprscan /home/linuxbrew/data/fun_newton/annotate_misc/iprscan.xml --busco_db nematoda



# 8) final annotate
funannotate annotate -i  ${out_folder} --cpus ${threads}


	






	
