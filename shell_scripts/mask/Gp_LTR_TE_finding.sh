#$ -cwd

cd /shelf/apps/User_name/newton/final_genome/

python run_pipeline.py


conda activate python27
conda activate repeatmasking

# genomtools is in path

/shelf/apps/User_name/apps/TransposonPSI_08222010/transposonPSI.pl Gp_Newton_haplotype1.fasta nucl

# LTRharvest in genome tools  http://genometools.org/documents/ltrharvest.pdf


#make the enhanced database
gt suffixerator -db Gp_Newton_haplotype1.fasta -indexname pea -tis -suf -lcp -des -ssp -sds -dna


#search for the beasties. 

gt ltrharvest -index pea -mintsd 5 -maxtsd 100 > pea_LTR_harvest.out


# (optional sequence clusterin gof output

#mkvtree .... look at the pdf. 
