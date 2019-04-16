
cd /storage/home/users/User_name/shelf_apps/newton/final_genome2


THREADS=12

conda activate pb_toolkit


conda activate haplo_phase
#minimap2 -t $THREADS -ax map-pb Gp_Newton_haplotype1.fasta \
#/shelf/apps/User_name/newton/reads/all_pacbio.fastq.gz > aln2.sam
#samtools view -@ $THREADS -S -b -o alnunsorted2.bam aln2.sam
#wait
#samtools sort -@ $THREADS -o  sorted_mini3.bam alnunsorted2.bam
#samtools index sorted_mini3.bam

#bwa index Gp_Newton_haplotype1.fasta
bwa mem -t $THREADS Gp_Newton_haplotype1.fasta \
/shelf/apps/User_name/newton/DNAseq/R1_prinseq_good_pm0d.fastq \
/shelf/apps/User_name/newton/DNAseq/R2_prinseq_good_vVZt.fastq > new.illumina.mapped2.sam

samtools view -@ $THREADS -S -b -o new.illumina.temp.mapped2.bam new.illumina.mapped2.sam
samtools sort -@ $THREADS -o  new.illumina.mapped2.bam new.illumina.temp.mapped2.bam
samtools index new.illumina.mapped2.bam


