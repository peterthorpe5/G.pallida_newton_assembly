


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

#funannotate mask -i Gpal_Lind_nanopore.fasta --cpus 40 -o Gpal_Lind_nanopore.mask.fasta --repeatmodeler_lib consensi.fa.censor

funannotate train --no_normalize_reads --trinity Trinity-GG.fasta --right_norm right.norm.fa --left_norm left.norm.fa -i Gpal_Lind_nanopore.mask.fasta   -o lindley --left R1.fq.gz --right R2.fq.gz --no_trimmomatic --memory 200G --species "gpal_lindly"  --cpus 40 --max_intronlen 10000
    
    
funannotate predict -i Gpal_Lind_nanopore.mask.fasta --name GPALLIND -o lindley --species  gpal_lindly --cpus 40 --max_intronlen 10000 --organism other --busco_db nematoda --database ../databases/ --organism other 

