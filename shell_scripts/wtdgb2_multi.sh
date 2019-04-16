#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
#$ -N MAIN_wtdbg ## job name
#$ -S /bin/bash ## shell where it will run this job
####$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory
###$ -q marvin.q ## queue name
#$ -pe multi 8
#$ -m e
#$ -M User_name@st-andrews.ac.uk

cd /storage/home/users/User_name/shelf_apps/newton

module load oraclejava/jdk1.8.0_74
#default

#/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -i Gp_newton_reduce_haplotypes_400X.correctedReads.fasta.gz -fo prefix

#/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i prefix.ctg.lay -fo prefix.ctg.lay.fa

#/shelf/apps/User_name/canu-1.7.1/Linux-amd64/bin/canu correctedErrorRate=0.85 corMhapSensitivity=normal -p newton_using_EC_reads_EC85 -d newton_using_EC_reads_EC85 genomeSize=120m -pacbio-raw 'Gp_newton_reduce_haplotypes_400X.correctedReads.fasta.gz' -useGrid=False -maxMemory=450 -maxThreads=16


#correctedErrorRate=0.085 corMhapSensitivity=normal

#mv ~/Hsac_all_data_reduce_ploidy_300X_CMAPSEN_normal.trimmedReads.fasta.gz ./


#/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -i Hsac_all_data_reduce_ploidy_300X_CMAPSEN_normal.trimmedReads.fasta.gz -fo Hs_trimmed

#/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i Hs_trimmed.ctg.lay -fo Hs_trimmed.ctg.lay.fa

#mv ~/*fastq.gz ./

#/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -i newton_all_PacBio.fastq.gz -fo newton_raw

#/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i newton_raw.ctg.lay -fo newton_raw.ctg.lay.fa

#/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -i all_pacbio.fastq.gz -fo Hs_raw

#/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i Hs_raw.ctg.lay -fo Hs_raw.ctg.lay.fa


#/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 --edge-min 2 --rescue-low-cov-edges -p 19 -i newton_using_EC_reads_EC85.correctedReads.fasta.gz -fo newt_correc_corre

#/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i newt_correc_corre.ctg.lay -fo newt_correc_corre.ctg.p19_rescuelay.fa

#/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 --edge-min 2 --rescue-low-cov-edges -p 19 -i Hsac_all_data_reduce_ploidy_300X_CMAPSEN_normal.trimmedReads.fasta.gz -fo Hs_trimmed

#/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i Hs_trimmed.ctg.lay -fo Hs_trimmed.ctg.p19_rescue.lay.fa


#mv ~/Hsac_all_data_default.correctedReads.fasta.gz ./


#/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 --edge-min 2 --rescue-low-cov-edges -p 19 -i Hsac_all_data_default.correctedReads.fasta.gz -fo Hs_default

#/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i Hs_default.ctg.lay -fo Hs_default.ctg.p19_rescue.lay.fa


#/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -i Hsac_all_data_default.correctedReads.fasta.gz -fo Hs_default1

#/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i Hs_default1.ctg.lay -fo Hs_default1.ctg.lay.fa


#/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -p 19 -L 5000 -i newton_using_EC_reads_EC85.correctedReads.fasta.
#gz -fo newt_correc_corre

#/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i newt_correc_corre.ctg.lay -fo newt_correc_corre.ctg.l5000_lay.fa


/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -p 19 -L 4000 -i newton_using_EC_reads_EC85.correctedReads.fasta.gz -fo newt_correc_corre4000

/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i newt_correc_corre4000.ctg.lay -fo newt_correc_corre.ctg.L4000_p19_lay.fa

/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -p 19 -L 5000 -i newton_using_EC_reads_EC85.correctedReads.fasta.gz -fo newt_correc_corre5000

/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i newt_correc_corre5000.ctg.lay -fo newt_correc_corre.ctg.L5000_p19_lay.fa

/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -L 5000 -i newton_using_EC_reads_EC85.correctedReads.fasta.gz -fo newt_correc_corre5000

/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i newt_correc_corre5000.ctg.lay -fo newt_correc_corre.ctg.L5000_lay.fa

/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -L 5000 -i Hsac_all_data_default.correctedReads.fasta.gz -fo Hs_default1

/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i Hs_default1.ctg.lay -fo Hs_default1_l5000.ctg.lay.fa



/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -L 5000 -p 19 -i Hsac_all_data_reduce_ploidy_300X_CMAPSEN_normal.trimmedReads.fasta.gz -fo Hs_trimmed

/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i Hs_trimmed.ctg.lay -fo Hs_trimmed.ctg.p19_L5000.lay.fa

/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -L 5000 -i newton_using_EC_reads_EC85.correctedReads.fasta.gz -fo newt_correc_corre1

/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i newt_correc_corre1.ctg.lay -fo newt_correc1_corre.ctg.L5000_p_default1_lay.fa

/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -L 5000 -i Hsac_all_data_reduce_ploidy_300X_CMAPSEN_normal.trimmedReads.fasta.gz -fo Hs_trimmed

/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i Hs_trimmed.ctg.lay -fo Hs_trimmed.ctg._L5000.lay.fa


/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -L 5000 -i Hsac_all_data_default.correctedReads.fasta.gz -fo Hs_trimmed

/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i Hs_trimmed.ctg.lay -fo Hs_default_trimmed.ctg._L5000.lay.fa

/shelf/apps/User_name/apps/wtdbg2-master/wtdbg2 -t 8 -L 5000 -i all_pacbio.fastq.gz -fo Hs_raw

/shelf/apps/User_name/apps/wtdbg2-master/wtpoa-cns -t 8 -i Hs_raw.ctg.lay -fo Hs_raw.ctg._L5000.lay.fa











