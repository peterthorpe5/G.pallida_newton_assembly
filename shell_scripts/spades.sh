cd /storage/home/users/User_name/shelf_apps/newton/DNAseq


module load SPAdes/3.12.0


spades.py -o Gp_newt -1 R1.fq.gz -2 R2.fq.gz -t 10 -m 250 

