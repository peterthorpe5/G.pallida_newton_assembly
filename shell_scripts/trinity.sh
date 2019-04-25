cd /storage/home/users/pjt6/shelf_apps/newton/RNAseq

conda activate trinity


 Trinity --genome_guided_bam GpAligned.sortedByCoord.out.bam \
         --genome_guided_max_intron 15000 \
         --max_memory 10G --CPU 12 --max_memory 100G --full_cleanup --output Gp_genome_gui_Trinity --genome_guided_min_coverage 5
         
 