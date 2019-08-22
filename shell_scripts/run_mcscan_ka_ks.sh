cd /storage/home/users/pjt6/synteny

cd hg_li
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/hg_li/hg_li 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g hg_li.gff -b hg_li.blast -c hg_li.collinearity -o hg_li_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g hg_li.gff -b hg_li.blast -c hg_li.collinearity -o hg_li_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i hg_li.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o hg_li.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd hg_gr
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/hg_gr/hg_gr 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g hg_gr.gff -b hg_gr.blast -c hg_gr.collinearity -o hg_gr_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g hg_gr.gff -b hg_gr.blast -c hg_gr.collinearity -o hg_gr_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i hg_gr.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o hg_gr.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd hg_mi
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/hg_mi/hg_mi 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g hg_mi.gff -b hg_mi.blast -c hg_mi.collinearity -o hg_mi_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g hg_mi.gff -b hg_mi.blast -c hg_mi.collinearity -o hg_mi_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i hg_mi.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o hg_mi.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd hg_gp
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/hg_gp/hg_gp 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g hg_gp.gff -b hg_gp.blast -c hg_gp.collinearity -o hg_gp_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g hg_gp.gff -b hg_gp.blast -c hg_gp.collinearity -o hg_gp_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i hg_gp.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o hg_gp.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd li_hg
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/li_hg/li_hg 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g li_hg.gff -b li_hg.blast -c li_hg.collinearity -o li_hg_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g li_hg.gff -b li_hg.blast -c li_hg.collinearity -o li_hg_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i li_hg.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o li_hg.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd li_gr
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/li_gr/li_gr 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g li_gr.gff -b li_gr.blast -c li_gr.collinearity -o li_gr_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g li_gr.gff -b li_gr.blast -c li_gr.collinearity -o li_gr_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i li_gr.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o li_gr.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd li_mi
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/li_mi/li_mi 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g li_mi.gff -b li_mi.blast -c li_mi.collinearity -o li_mi_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g li_mi.gff -b li_mi.blast -c li_mi.collinearity -o li_mi_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i li_mi.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o li_mi.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd li_gp
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/li_gp/li_gp 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g li_gp.gff -b li_gp.blast -c li_gp.collinearity -o li_gp_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g li_gp.gff -b li_gp.blast -c li_gp.collinearity -o li_gp_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i li_gp.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o li_gp.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd gr_hg
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/gr_hg/gr_hg 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gr_hg.gff -b gr_hg.blast -c gr_hg.collinearity -o gr_hg_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gr_hg.gff -b gr_hg.blast -c gr_hg.collinearity -o gr_hg_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i gr_hg.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o gr_hg.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd gr_li
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/gr_li/gr_li 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gr_li.gff -b gr_li.blast -c gr_li.collinearity -o gr_li_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gr_li.gff -b gr_li.blast -c gr_li.collinearity -o gr_li_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i gr_li.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o gr_li.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd gr_mi
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/gr_mi/gr_mi 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gr_mi.gff -b gr_mi.blast -c gr_mi.collinearity -o gr_mi_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gr_mi.gff -b gr_mi.blast -c gr_mi.collinearity -o gr_mi_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i gr_mi.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o gr_mi.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd gr_gp
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/gr_gp/gr_gp 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gr_gp.gff -b gr_gp.blast -c gr_gp.collinearity -o gr_gp_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gr_gp.gff -b gr_gp.blast -c gr_gp.collinearity -o gr_gp_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i gr_gp.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o gr_gp.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd mi_hg
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/mi_hg/mi_hg 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g mi_hg.gff -b mi_hg.blast -c mi_hg.collinearity -o mi_hg_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g mi_hg.gff -b mi_hg.blast -c mi_hg.collinearity -o mi_hg_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i mi_hg.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o mi_hg.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd mi_li
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/mi_li/mi_li 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g mi_li.gff -b mi_li.blast -c mi_li.collinearity -o mi_li_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g mi_li.gff -b mi_li.blast -c mi_li.collinearity -o mi_li_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i mi_li.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o mi_li.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd mi_gr
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/mi_gr/mi_gr 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g mi_gr.gff -b mi_gr.blast -c mi_gr.collinearity -o mi_gr_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g mi_gr.gff -b mi_gr.blast -c mi_gr.collinearity -o mi_gr_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i mi_gr.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o mi_gr.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd mi_gp
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/mi_gp/mi_gp 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g mi_gp.gff -b mi_gp.blast -c mi_gp.collinearity -o mi_gp_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g mi_gp.gff -b mi_gp.blast -c mi_gp.collinearity -o mi_gp_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i mi_gp.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o mi_gp.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd gp_hg
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/gp_hg/gp_hg 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gp_hg.gff -b gp_hg.blast -c gp_hg.collinearity -o gp_hg_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gp_hg.gff -b gp_hg.blast -c gp_hg.collinearity -o gp_hg_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i gp_hg.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o gp_hg.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd gp_li
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/gp_li/gp_li 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gp_li.gff -b gp_li.blast -c gp_li.collinearity -o gp_li_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gp_li.gff -b gp_li.blast -c gp_li.collinearity -o gp_li_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i gp_li.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o gp_li.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd gp_gr
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/gp_gr/gp_gr 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gp_gr.gff -b gp_gr.blast -c gp_gr.collinearity -o gp_gr_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gp_gr.gff -b gp_gr.blast -c gp_gr.collinearity -o gp_gr_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i gp_gr.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o gp_gr.collinearity.kaks 
cd /storage/home/users/pjt6/synteny
cd gp_mi
~/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/pjt6/synteny/gp_mi/gp_mi 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gp_mi.gff -b gp_mi.blast -c gp_mi.collinearity -o gp_mi_detect_collinear_tandem_array 
~/shelf_apps/apps/MCScanX/downstream_analyses/detect_collinear_tandem_arrays -g gp_mi.gff -b gp_mi.blast -c gp_mi.collinearity -o gp_mi_detect_collinear_tandem_array.out 
perl ~/shelf_apps/apps/MCScanX/downstream_analyses/add_ka_and_ks_to_collinearity.pl -i gp_mi.collinearity -d /storage/home/users/pjt6/synteny/all_nt.fasta -o gp_mi.collinearity.kaks 
cd /storage/home/users/pjt6/synteny