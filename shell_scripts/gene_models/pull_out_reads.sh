cd /storage/home/users/User_name/newton/gene_models/brakerrnaseq/transrate_results_done/nt

module load samtools
samtools view -b -f 4 R1.fq.gz.R2.fq.gz.nt.bam > unmapped.bam



#conda activate bedtools

samtools sort -@ 4 -n unmapped.bam aln.qsort


bedtools bamtofastq -i aln.qsort.bam \
                      -fq unmapped_R1.fq \
                      -fq2 unmapped_R2.fq
