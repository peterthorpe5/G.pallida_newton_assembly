cd /storage/home/users/pjt6/newton/spry/w_Hsch/spyrsecs/fasta_files/



# no trim al
/storage/home/users/pjt6/shelf_apps/apps/iqtree-1.6.12-Linux/bin/iqtree -s all_w_hs.SPRYCES_min70.hmm.fasta \
-st AA -m TEST -bb 1000 -alrt 1000 -mem 1250GB -nt AUTO 

# trim al gappy alignment
/storage/home/users/pjt6/shelf_apps/apps/iqtree-1.6.12-Linux/bin/iqtree -s all_w_hs.SPRYCES_min70.hmm.gappy.trimal.fasta \
-st AA -m TEST -bb 1000 -alrt 1000 -mem 1250GB -nt AUTO 

