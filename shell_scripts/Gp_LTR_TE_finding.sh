#$ -cwd

cd /shelf/apps/User_name/newton/final_genome2/

#python run_pipeline.py


conda activate python27
conda activate repeatmasking

# genomtools is in path

#/shelf/apps/User_name/apps/TransposonPSI_08222010/transposonPSI.pl Gp_Newton_haplotype1.fasta nuc

# LTRharvest in genome tools  http://genometools.org/documents/ltrharvest.pdf


#make the enhanced database
gt suffixerator -db Gp_Newton_haplotype1.fasta -indexname gp -tis -suf -lcp -des -ssp -sds -dna


#search for the beasties. 

gt ltrharvest -index gp -mintsd 5 -maxtsd 100 > Gp_LTR_harvest.out


# (optional sequence clusterin gof output

#mkvtree .... look at the pdf. 
