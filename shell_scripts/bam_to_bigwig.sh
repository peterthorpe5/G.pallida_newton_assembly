cd /storage/home/users/pjt6/project/Lindley_RNAseq_mapping
	
bedtools genomecov -bga -ibam Gp_14DPI_rep1Aligned.sortedByCoord.out.bam  > Gp_14DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_14DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph >Gp_14DPI_rep1Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_14DPI_rep1Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_14DPI_rep1Aligned.sortedByCoord.out.bam.bw   && rm Gp_14DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_7DPI_rep2Aligned.sortedByCoord.out.bam  > Gp_7DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_7DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph >Gp_7DPI_rep2Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_7DPI_rep2Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_7DPI_rep2Aligned.sortedByCoord.out.bam.bw   && rm Gp_7DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_14DPI_rep2Aligned.sortedByCoord.out.bam  > Gp_14DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_14DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph >Gp_14DPI_rep2Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_14DPI_rep2Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_14DPI_rep2Aligned.sortedByCoord.out.bam.bw   && rm Gp_14DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_EGG_rep1Aligned.sortedByCoord.out.bam  > Gp_EGG_rep1Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_EGG_rep1Aligned.sortedByCoord.out.bam.bedgraph >Gp_EGG_rep1Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_EGG_rep1Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_EGG_rep1Aligned.sortedByCoord.out.bam.bw   && rm Gp_EGG_rep1Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_21DPI_rep1Aligned.sortedByCoord.out.bam  > Gp_21DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_21DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph >Gp_21DPI_rep1Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_21DPI_rep1Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_21DPI_rep1Aligned.sortedByCoord.out.bam.bw   && rm Gp_21DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_EGG_rep2Aligned.sortedByCoord.out.bam  > Gp_EGG_rep2Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_EGG_rep2Aligned.sortedByCoord.out.bam.bedgraph >Gp_EGG_rep2Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_EGG_rep2Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_EGG_rep2Aligned.sortedByCoord.out.bam.bw   && rm Gp_EGG_rep2Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_21DPI_rep2Aligned.sortedByCoord.out.bam  > Gp_21DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_21DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph >Gp_21DPI_rep2Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_21DPI_rep2Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_21DPI_rep2Aligned.sortedByCoord.out.bam.bw   && rm Gp_21DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_J2_rep1Aligned.sortedByCoord.out.bam  > Gp_J2_rep1Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_J2_rep1Aligned.sortedByCoord.out.bam.bedgraph >Gp_J2_rep1Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_J2_rep1Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_J2_rep1Aligned.sortedByCoord.out.bam.bw   && rm Gp_J2_rep1Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_28DPI_rep1Aligned.sortedByCoord.out.bam  > Gp_28DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_28DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph >Gp_28DPI_rep1Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_28DPI_rep1Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_28DPI_rep1Aligned.sortedByCoord.out.bam.bw   && rm Gp_28DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_J2_rep2Aligned.sortedByCoord.out.bam  > Gp_J2_rep2Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_J2_rep2Aligned.sortedByCoord.out.bam.bedgraph >Gp_J2_rep2Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_J2_rep2Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_J2_rep2Aligned.sortedByCoord.out.bam.bw   && rm Gp_J2_rep2Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_28DPI_rep2Aligned.sortedByCoord.out.bam  > Gp_28DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_28DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph >Gp_28DPI_rep2Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_28DPI_rep2Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_28DPI_rep2Aligned.sortedByCoord.out.bam.bw   && rm Gp_28DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_J2_rep3Aligned.sortedByCoord.out.bam  > Gp_J2_rep3Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_J2_rep3Aligned.sortedByCoord.out.bam.bedgraph >Gp_J2_rep3Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_J2_rep3Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_J2_rep3Aligned.sortedByCoord.out.bam.bw   && rm Gp_J2_rep3Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_35DPI_rep1Aligned.sortedByCoord.out.bam  > Gp_35DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_35DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph >Gp_35DPI_rep1Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_35DPI_rep1Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_35DPI_rep1Aligned.sortedByCoord.out.bam.bw   && rm Gp_35DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_MALE_rep1_cAligned.sortedByCoord.out.bam  > Gp_MALE_rep1_cAligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_MALE_rep1_cAligned.sortedByCoord.out.bam.bedgraph >Gp_MALE_rep1_cAligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_MALE_rep1_cAligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_MALE_rep1_cAligned.sortedByCoord.out.bam.bw   && rm Gp_MALE_rep1_cAligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_35DPI_rep2Aligned.sortedByCoord.out.bam  > Gp_35DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_35DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph >Gp_35DPI_rep2Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_35DPI_rep2Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_35DPI_rep2Aligned.sortedByCoord.out.bam.bw   && rm Gp_35DPI_rep2Aligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_MALE_rep2_bAligned.sortedByCoord.out.bam  > Gp_MALE_rep2_bAligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_MALE_rep2_bAligned.sortedByCoord.out.bam.bedgraph >Gp_MALE_rep2_bAligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_MALE_rep2_bAligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_MALE_rep2_bAligned.sortedByCoord.out.bam.bw   && rm Gp_MALE_rep2_bAligned.sortedByCoord.out.bam.bedgraph 

bedtools genomecov -bga -ibam Gp_7DPI_rep1Aligned.sortedByCoord.out.bam  > Gp_7DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph  && bedtools sort -i  Gp_7DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph >Gp_7DPI_rep1Aligned.sortedByCoord.out.bam2.bedgraph && bedGraphToBigWig Gp_7DPI_rep1Aligned.sortedByCoord.out.bam2.bedgraph test.fasta.chrom.sizes Gp_7DPI_rep1Aligned.sortedByCoord.out.bam.bw   && rm Gp_7DPI_rep1Aligned.sortedByCoord.out.bam.bedgraph 

samtools index merged.bam
  bedtools genomecov -bga -ibam merged.bam  > merged.bam.bedgraph  && bedtools sort -i  merged.bam.bedgraph >merged.bam2.bedgraph && bedGraphToBigWig merged.bam2.bedgraph test.fasta.chrom.sizes merged.bam.bw   && rm merged.bam.bedgraph

# samtols merge no header

 samtools merge -@ 8 merged_no_head  Gp_14DPI_rep1Aligned.sortedByCoord.out.bam  Gp_14DPI_rep2Aligned.sortedByCoord.out.bam  Gp_21DPI_rep1Aligned.sortedByCoord.out.bam  Gp_21DPI_rep2Aligned.sortedByCoord.out.bam  Gp_28DPI_rep1Aligned.sortedByCoord.out.bam  Gp_28DPI_rep2Aligned.sortedByCoord.out.bam  Gp_35DPI_rep1Aligned.sortedByCoord.out.bam  Gp_35DPI_rep2Aligned.sortedByCoord.out.bam  Gp_7DPI_rep1Aligned.sortedByCoord.out.bam  Gp_7DPI_rep2Aligned.sortedByCoord.out.bam  Gp_EGG_rep1Aligned.sortedByCoord.out.bam  Gp_EGG_rep2Aligned.sortedByCoord.out.bam  Gp_J2_rep1Aligned.sortedByCoord.out.bam  Gp_J2_rep2Aligned.sortedByCoord.out.bam  Gp_J2_rep3Aligned.sortedByCoord.out.bam  Gp_MALE_rep1_cAligned.sortedByCoord.out.bam  Gp_MALE_rep2_bAligned.sortedByCoord.out.bam
samtools index  merged_no_head

 bedtools genomecov -bga -imerged_no_head.bam Gp_21DPI_rep1Aligned.sortedByCoord.out.merged_no_head.bam  > Gp_21DPI_rep1Aligned.sortedByCoord.out.merged_no_head.bam.bedgraph  && bedtools sort -i  Gp_21DPI_rep1Aligned.sortedByCoord.out.merged_no_head.bam.bedgraph >Gp_21DPI_rep1Aligned.sortedByCoord.out.merged_no_head.bam2.bedgraph && bedGraphToBigWig Gp_21DPI_rep1Aligned.sortedByCoord.out.merged_no_head.bam2.bedgraph test.fasta.chrom.sizes Gp_21DPI_rep1Aligned.sortedByCoord.out.merged_no_head.bam.bw   && rm Gp_21DPI_rep1Aligned.sortedByCoord.out.merged_no_head.bam.bedgraph
