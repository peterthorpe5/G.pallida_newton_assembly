
cd /shelf/apps/User_name/newton/stricter_scaff_polish 

module load hmmer/3.1b2

hmmsearch --cut_ga --domtblout Gp_newton_gene_vs_spry_domain.out hmm gpal_lindley.proteins.fa


perl /storage/home/users/pjt6/shelf_apps/apps/phobius/phobius.pl -short gpal_lindley.proteins.fa > Gp_newton.phobius


cat Gp_newton_gene_vs_spry_domain.out | grep -v "#" |  tr -s ' ' | cut -d ' ' -f1 | grep "y1"> SPRY.names

cat  Gp_newton.phobius | tr -s ' ' | grep "0 Y" | cut -d ' ' -f1 > secreted.names

