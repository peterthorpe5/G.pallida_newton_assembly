cd /storage/home/users/User_name/shelf_apps/newton/final_genome2
bedtools genomecov -ibam sorted_mini3.bam -bga | awk '$4==0' > assembly_zero_cov.txt

