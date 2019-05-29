cd /storage/home/users/pjt6/newton/final_genes/funaannot/update_results/effector_expression_clustering

#conda activate trinity
#cd Gp_putative_effectors

 cat Gp_genes.counts.matrix.TMM_normalized.FPKM | head -n 1 > Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM

 cat Gp_genes.counts.matrix.TMM_normalized.FPKM | grep -Ff all_putative_effectors.txt >> Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM

# coloumn oderering
/shelf/apps/pjt6/apps/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/PtR --order_columns_by_samples_file --heatmap_scale_limits "-5,5" -s samples_described.txt -m Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM --save --heatmap --log2 --center_rows



/shelf/apps/pjt6/apps/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --no_column_reordering --Ptree 10 -R Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM.RData
mv clusters_fixed_* ./Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM.RData.clusters_fixed_P_10

/shelf/apps/pjt6/apps/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 20 -R Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM.RData
mv clusters_fixed_* ./Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM.RData.clusters_fixed_P_20
/shelf/apps/pjt6/apps/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 30 -R Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM.RData
mv clusters_fixed_* ./Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM.RData.clusters_fixed_P_30
/shelf/apps/pjt6/apps/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM.RData
mv clusters_fixed_* ./Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM.RData.clusters_fixed_P_50
/shelf/apps/pjt6/apps/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering  --Ptree 70 -R Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM.RData
mv clusters_fixed_* ./Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM.RData.clusters_fixed_P_70

/shelf/apps/pjt6/apps/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering  --Ktree 8 -R Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM.RData
mv clusters_fixed_Ktree_* ./Gp_putative_effectors.genes.counts.matrix.TMM_normalized.FPKM.RData.clusters_fixed_Ktree_8

