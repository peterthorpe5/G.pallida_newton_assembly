cd /home/User_name/Desktop/newton_final


#generate protein alingmnert:
#gth -genomic Gp_Newton_haplotype1_sotmasked.fasta -protein Gpal_aa_fromNCBI_JJ.fasta -skipalignmentout -gff3out -o gth.aln -force -paralogs -prseedlength 20 -prhdist 2 -gcmincoverage 80 -prminmatchlen 20

#############################################################
# 1) just braker nothing extra
cd /home/User_name/Desktop/newton_final/braker_20_rounds

mv ../Gp_Newton_haplotype1_sotmasked.fasta ./
mv ../Gp_soft_maskedAligned.sortedByCoord.out.bam ./

echo "1) running test on comp braker and utr nothing extra"
cmd='perl /home/User_name/programs/BRAKER/scripts/braker.pl 
--genome=Gp_Newton_haplotype1_sotmasked.fasta
--overwrite 
--workingdir=/home/User_name/Desktop/newton_final/braker_20_rounds 
--bam=Gp_soft_maskedAligned.sortedByCoord.out.bam 
--AUGUSTUS_BIN_PATH=/home/User_name/programs/Augustus/bin 
--AUGUSTUS_CONFIG_PATH=/home/User_name/programs/Augustus/config 
--AUGUSTUS_SCRIPTS_PATH=/home/User_name/programs/Augustus/scripts 
--species=GPAL_newton_softmasked20rounds0
--augustus_args="--protein=on --start=on --stop=on --cds=on --introns=on  --stopCodonExcludedFromCDS=False --genemodel=complete " 
--crf --cores 4 --gff3 --rounds 20 
--filterOutShort 
--GENEMARK_PATH=/home/User_name/programs/gm_et_linux_64/gmes_petap/ 
--softmasking '
echo ${cmd}
eval ${cmd}

mv ./Gp_Newton_haplotype1_sotmasked.fasta ../
mv ./Gp_soft_maskedAligned.sortedByCoord.out.bam ../



#############################################################
# 1b) just braker nothing extra - with stop
cd /home/User_name/Desktop/newton_final/braker_stop 

mv ../Gp_Newton_haplotype1_sotmasked.fasta ./
mv ../Gp_soft_maskedAligned.sortedByCoord.out.bam ./

echo "1) running test on comp braker and utr nothing extra"
cmd='perl /home/User_name/programs/BRAKER/scripts/braker.pl 
--genome=Gp_Newton_haplotype1_sotmasked.fasta
--overwrite 
--workingdir=/home/User_name/Desktop/newton_final/braker 
--bam=Gp_soft_maskedAligned.sortedByCoord.out.bam 
--AUGUSTUS_BIN_PATH=/home/User_name/programs/Augustus/bin 
--AUGUSTUS_CONFIG_PATH=/home/User_name/programs/Augustus/config 
--AUGUSTUS_SCRIPTS_PATH=/home/User_name/programs/Augustus/scripts 
--species=GPAL_newton_softmasked20stop
--augustus_args=" --stopCodonExcludedFromCDS=False --genemodel=complete  " 
--crf --cores 4 --gff3 --rounds 20 
--utr on
--filterOutShort 
--GENEMARK_PATH=/home/User_name/programs/gm_et_linux_64/gmes_petap/ 
--softmasking '
echo ${cmd}
eval ${cmd}

mv ./Gp_Newton_haplotype1_sotmasked.fasta ../
mv ./Gp_soft_maskedAligned.sortedByCoord.out.bam ../



#############################################################
# 2)  braker and utr predixction
cd /home/User_name/Desktop/newton_final/utr

mv ../Gp_Newton_haplotype1_sotmasked.fasta ./
mv ../Gp_soft_maskedAligned.sortedByCoord.out.bam ./

echo "running test on comp braker and utr nothing extra"
cmd='perl /home/User_name/programs/BRAKER/scripts/braker.pl 
--genome=Gp_Newton_haplotype1_sotmasked.fasta
--overwrite 
--workingdir=/home/User_name/Desktop/newton_final/utr 
--bam=Gp_soft_maskedAligned.sortedByCoord.out.bam 
--AUGUSTUS_BIN_PATH=/home/User_name/programs/Augustus/bin 
--AUGUSTUS_CONFIG_PATH=/home/User_name/programs/Augustus/config 
--AUGUSTUS_SCRIPTS_PATH=/home/User_name/programs/Augustus/scripts 
--species=testing_UTR30_w_cds
--augustus_args="--protein=on --start=on --stop=on --cds=on --introns=on --noInFrameStop=true --genemodel=complete " 
--crf --cores 4 --gff3 --rounds 10 
--filterOutShort
--UTR=on
--GENEMARK_PATH=/home/User_name/programs/gm_et_linux_64/gmes_petap/ 
--softmasking '
#echo ${cmd}
#eval ${cmd}

mv ./Gp_Newton_haplotype1_sotmasked.fasta ../
mv ./Gp_soft_maskedAligned.sortedByCoord.out.bam ../


#############################################################
# 2b)  braker and utr predixction with stop for funannotate
cd /home/User_name/Desktop/newton_final/utr_stop

mv ../Gp_Newton_haplotype1_sotmasked.fasta ./
mv ../Gp_soft_maskedAligned.sortedByCoord.out.bam ./

echo "running test on comp braker and utr nothing extra"
cmd='perl /home/User_name/programs/BRAKER/scripts/braker.pl 
--genome=Gp_Newton_haplotype1_sotmasked.fasta
--overwrite 
--workingdir=/home/User_name/Desktop/newton_final/utr 
--bam=Gp_soft_maskedAligned.sortedByCoord.out.bam 
--AUGUSTUS_BIN_PATH=/home/User_name/programs/Augustus/bin 
--AUGUSTUS_CONFIG_PATH=/home/User_name/programs/Augustus/config 
--AUGUSTUS_SCRIPTS_PATH=/home/User_name/programs/Augustus/scripts 
--species=testing_UTR30_w_cds
--augustus_args=" --stopCodonExcludedFromCDS=False --genemodel=complete  " 
--crf --cores 4 --gff3 --rounds 10 
--filterOutShort
--UTR=on
--GENEMARK_PATH=/home/User_name/programs/gm_et_linux_64/gmes_petap/ 
--softmasking '
echo ${cmd}
eval ${cmd}

mv ./Gp_Newton_haplotype1_sotmasked.fasta ../
mv ./Gp_soft_maskedAligned.sortedByCoord.out.bam ../



#############################################################
# 3)  braker and utr predixction with proteins
cd /home/User_name/Desktop/newton_final/with_UTR_cds_info


mv ../Gp_Newton_haplotype1_sotmasked.fasta ./
mv ../Gp_soft_maskedAligned.sortedByCoord.out.bam ./


echo "running test on comp - predict utr WITH POTEIN AS INPUT additional trainin with proteins"
cmd='perl /home/User_name/programs/BRAKER/scripts/braker.pl 
--genome=Gp_Newton_haplotype1_sotmasked.fasta
--overwrite 
--workingdir=/home/User_name/Desktop/newton_final/with_UTR_cds_info 
--bam=Gp_soft_maskedAligned.sortedByCoord.out.bam 
--AUGUSTUS_BIN_PATH=/home/User_name/programs/Augustus/bin 
--AUGUSTUS_CONFIG_PATH=/home/User_name/programs/Augustus/config 
--AUGUSTUS_SCRIPTS_PATH=/home/User_name/programs/Augustus/scripts 
--species=testing_UTR_w_cds
--gth2traingenes
--prot_aln=gth.aln
--prg=gth
--ALIGNMENT_TOOL_PATH=~/programs/gth-1.7.1-Linux_x86_64-64bit/bin/
--augustus_args="--protein=on --start=on --stop=on --cds=on --introns=on --noInFrameStop=true --genemodel=complete " 
--crf --cores 4 --gff3 --rounds 1 
--UTR=on --filterOutShort 
--GENEMARK_PATH=/home/User_name/programs/gm_et_linux_64/gmes_petap/ 
--softmasking '
echo ${cmd}
eval ${cmd}


mv ./Gp_Newton_haplotype1_sotmasked.fasta ../
mv ./Gp_soft_maskedAligned.sortedByCoord.out.bam ../


#############################################################
# 4)  braker NO utr predixction with proteins

cd /home/User_name/Desktop/newton_final/prot_no_utr

mv ../Gp_Newton_haplotype1_sotmasked.fasta ./
mv ../Gp_soft_maskedAligned.sortedByCoord.out.bam ./

echo "running test on comp braker and prot nothing extra"
cmd='perl /home/User_name/programs/BRAKER/scripts/braker.pl 
--genome=Gp_Newton_haplotype1_sotmasked.fasta
--overwrite 
--workingdir=/home/User_name/Desktop/newton_final/prot_no_utr 
--bam=Gp_soft_maskedAligned.sortedByCoord.out.bam 
--AUGUSTUS_BIN_PATH=/home/User_name/programs/Augustus/bin 
--AUGUSTUS_CONFIG_PATH=/home/User_name/programs/Augustus/config 
--AUGUSTUS_SCRIPTS_PATH=/home/User_name/programs/Augustus/scripts 
--species=testing_UTR40_w_cds
--augustus_args="--protein=on --start=on --stop=on --cds=on --introns=on --noInFrameStop=true --genemodel=complete " 
--crf --cores 4 --gff3 --rounds 10 
--filterOutShort
--gth2traingenes
--prot_aln=gth.aln
--prg=gth
--ALIGNMENT_TOOL_PATH=~/programs/gth-1.7.1-Linux_x86_64-64bit/bin/
--GENEMARK_PATH=/home/User_name/programs/gm_et_linux_64/gmes_petap/ 
--softmasking '
echo ${cmd}
eval ${cmd}


mv ./Gp_Newton_haplotype1_sotmasked.fasta ../
mv ./Gp_soft_maskedAligned.sortedByCoord.out.bam ../

#############################################################
# 4b)  braker NO utr predixction with proteins

cd /home/User_name/Desktop/newton_final/prot_no_utr_stop

mv ../Gp_Newton_haplotype1_sotmasked.fasta ./
mv ../Gp_soft_maskedAligned.sortedByCoord.out.bam ./

echo "running test on comp braker and prot nothing extra"
cmd='perl /home/User_name/programs/BRAKER/scripts/braker.pl 
--genome=Gp_Newton_haplotype1_sotmasked.fasta
--overwrite 
--workingdir=/home/User_name/Desktop/newton_final/prot_no_utr 
--bam=Gp_soft_maskedAligned.sortedByCoord.out.bam 
--AUGUSTUS_BIN_PATH=/home/User_name/programs/Augustus/bin 
--AUGUSTUS_CONFIG_PATH=/home/User_name/programs/Augustus/config 
--AUGUSTUS_SCRIPTS_PATH=/home/User_name/programs/Augustus/scripts 
--species=testing_UTR40_w_cds
--augustus_args="--protein=on --start=on --stop=on --cds=on --introns=on --noInFrameStop=true --genemodel=complete " 
--crf --cores 4 --gff3 --rounds 10 
--filterOutShort
--gth2traingenes
--prot_aln=gth.aln
--prg=gth
--ALIGNMENT_TOOL_PATH=~/programs/gth-1.7.1-Linux_x86_64-64bit/bin/
--GENEMARK_PATH=/home/User_name/programs/gm_et_linux_64/gmes_petap/ 
--softmasking '
echo ${cmd}
eval ${cmd}

mv ./Gp_Newton_haplotype1_sotmasked.fasta ../
mv ./Gp_soft_maskedAligned.sortedByCoord.out.bam ../

####################################################################################################
###################################################################################################
# gff foratm for funanaotate


# " --stopCodonExcludedFromCDS=False" command for funnatate augustus gff

