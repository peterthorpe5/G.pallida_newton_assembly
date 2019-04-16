#!/bin/bash
#$ -cwd
cd ~/shelf_apps/newton/stricter_scaff_polish

conda activate python27


minialign -t16 -xont Newton_Illumina_pilon_iter5.fasta ../newton_all_PacBio.fastq.gz > temp.sam
wait 
#samtools view -@ 18 -S -b -o unsorted.bam temp.sam
wait 
#samtools sort -@ 32 unsorted.bam sorted

#samtools index sorted

wait 
#purge_haplotigs readhist sorted.bam

#purge_haplotigs  contigcov  -i sorted.bam.genecov -o coverage_stats.csv  -l 3  -m 27  -h 130

#

#purge_haplotigs purge  -g Newton_Illumina_pilon_iter5.fasta  -c coverage_stats.csv  -b sorted.bam  -t 18  -a 60

#cp curated.haplotigs.fasta reassigned.haplotigs.fasta


#python ~/Downloads/blobtools-light-master/mapping2cov.py -a /home/pt40963/scratch/nematode/newton/purge_haplo_ALL_haplotypes_all_data/purge_the_purged/Newton_Illumina_pilon_iter5.fasta -bam /home/pt40963/scratch/nematode/newton/purge_haplo_ALL_haplotypes_all_data/purge_the_purged/sorted.bam -o contig.mappingall_final.cas.cov


######################################
blastn -db globodera_rostochiensis.PRJEB13504.WBPS12.genomic.fa -query Newton_Illumina_pilon_iter5.fasta -outfmt \
'6 qseqid staxids bitscore std scomnames sscinames sblastnames sskingdoms stitle' \
-evalue 1e-30 -out n.vGROS.out -num_threads 16

blastn -task megablast -query Newton_Illumina_pilon_iter5.fasta -db nt -outfmt '6 qseqid staxids bitscore std scomnames sscinames sblastnames sskingdoms stitle' \
-evalue 1e-30 -out n.clc.allfinal.out -num_threads 16

blobtools create -i Newton_Illumina_pilon_iter5.fasta  -s temp.sam -t n.clc.allfinal.out -t n.vGROS.out -o Newton_Illumina_pilon_iter5.fasta.blobplots
 
 mkdir allfinal.fa.blobplots
 cp Newton_Illumina_pilon_iter5.fasta.blobplots.blobDB.json ./allfinal.fa.blobplots
 
blobtools view -i  Newton_Illumina_pilon_iter5.fasta.blobplots.blobDB.json -o allfinal.fa.blobplots/
 
blobtools blobplot -i Newton_Illumina_pilon_iter5.fasta.blobplots.blobDB.json -o allfinal.fa.blobplots/ 

#rm temp.sam

cd allfinal.fa.blobplots
cat *.blobDB.table.txt | grep "Streptophyta" | cut -f1 > Streptophyta.names
#cat *.blobDB.table.txt | grep "Arthropoda" | cut -f1 > Arthropoda.names
cat *.blobDB.table.txt | grep "Ascomycota" | cut -f1 > Ascomycota.names
#cat *.blobDB.table.txt | grep "Chordata" | cut -f1 > Chordata.names
#cat *.blobDB.table.txt | grep "Nematoda" | cut -f1 > Nematoda.names
cat *.blobDB.table.txt | grep "Basidiomycota" | cut -f1 > Basidiomycota.names
cat *.blobDB.table.txt | grep "Proteobacteria" | cut -f1 > Proteobacteria.names
cat *.blobDB.table.txt | grep "Bacteria-undef" | cut -f1 > Bacteria-undef.names
cat *.blobDB.table.txt | grep "Viruses-undef" | cut -f1 > Viruses-undef.names
#cat *.blobDB.table.txt | grep "Cnidaria" | cut -f1 > Cnidaria.names
cat *.blobDB.table.txt | grep "Actinobacteria" | cut -f1 > Actinobacteria.names
cat *.blobDB.table.txt | grep "Mucoromycota" | cut -f1 > Mucoromycota.names
cat *.blobDB.table.txt | grep "Euryarchaeota" | cut -f1 > Euryarchaeota.names 
cat *.blobDB.table.txt | awk '$5 < 5' | grep 'no-hit' | sort -k4 -rn | cut -f1 > low_cov_5.LOOKATME
cat *.names >  bad_contigs.out

cd ../

python /shelf/apps/User_name/apps/public_scripts-master/get_sequences_i_want_from_fasta_command_line_not_wanted_file.py Newton_Illumina_pilon_iter5.fasta ./allfinal.fa.blobplots/bad_contigs.out 10 Hs_default_trimmed.ctg._L5000.lay.contim_filtered.fasta
rm temp.sam

#cp newt_correc_corre.ctg.L4000_p19.contim_filtered_final.fasta contigs.fasta

#python /shelf/apps/User_name/apps/finishingTool/finisherSC.py -par 8 /storage/home/users/User_name/shelf_apps/newton/newton_wtdgb/L4000/ /shelf/apps/User_name/conda/envs/python27/bin/


