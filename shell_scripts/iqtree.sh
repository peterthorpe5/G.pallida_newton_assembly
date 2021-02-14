cd /storage/home/users/pjt6/newton/spry
# trim al gappy alignment
#/storage/home/users/pjt6/shelf_apps/apps/iqtree-1.6.12-Linux/bin/iqtree -s al_nems_vs_SPRY_domains_min_70AA_trimAl.fasta \
-st AA -m TEST -bb 1000 -alrt 1000 -mem 200GB -nt AUTO 

# no trim al
/storage/home/users/pjt6/shelf_apps/apps/iqtree-1.6.12-Linux/bin/iqtree -s /storage/home/users/pjt6/newton/spry/fasta_files/al_nems_vs_SPRY_domains_min_70AA.fastaaln.fasta \
-st AA -m TEST -bb 1000 -alrt 1000 -mem 200GB -nt AUTO 
