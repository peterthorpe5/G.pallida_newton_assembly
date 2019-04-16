
cd /shelf/apps/User_name/newton/stricter_scaff_polish 

#module load hmmer/3.1b2

#hmmsearch --cut_ga --domtblout Gp_newton_gene_vs_spry_domain.out hmm braker_predictions.aa


perl /storage/home/users/User_name/shelf_apps/apps/phobius/phobius.pl -short \
braker_predictions_no_stop.aa > Gp_newton.phobius