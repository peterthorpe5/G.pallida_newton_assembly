

cd /storage/home/users/User_name/medicago/A17/DE_analysis_kallisto


#Abort on any error,
set -e

#conda activate trinity

cuwodi=/storage/home/users/User_name/medicago/A17/DE_analysis_kallisto

outdir=DE_analysis_EdgeRLOG2

cd ${cuwodi}/


TRINITY_HOME=/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.8.4/



 
 


#echo "\tGp_7DPI_rep1\tGp_7DPI_rep2\tGp_14DPI_rep1\tGp_14DPI_rep2\tGp_21DPI_rep1\tGp_21DPI_rep2\tGp_28DPI_rep1\tGp_28DPI_rep2\tGp_35DPI_rep1\tGp_35DPI_rep2\tGp_EGG_rep1\tGp_EGG_rep2\tGp_J2_rep1\tGp_J2_rep2\tGp_J2_rep3\tGp_MALE_rep1\tGp_MALE_rep2" > Gp_genes.counts.matrix
# merge these into a file
#FILES=$(ls -t -v *.counts | tr '\n' ' ')
#awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES | grep -v "__" >> Gp_genes.counts.matrix


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix *genes.counts.matrix --samples_file samples_described.txt --method edgeR --output ${outdir}
wait 


wait 


#/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --matrix *genes.counts.matrix --lengths genes_feature_lengths.txt


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --lexical_column_ordering --matrix *genes.counts.matrix --samples samples_described.txt --CPM --log2 --min_rowSums 10 --compare_replicates

wait
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --lexical_column_ordering --matrix *genes.counts.matrix --min_rowSums 2 --samples samples_described.txt --log2 --CPM --sample_cor_matrix

wait
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix *genes.counts.matrix -s samples_described.txt --min_rowSums 10 --log2 --CPM --center_rows --prin_comp 3

# does not need gene lenghts for this method
#/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --matrix *genes.counts.matrix --just_do_TMM_scaling


# .genes.counts.matrix.TMM_normalized.FPKM_normalized.FPKM - if you use gene lengths
# *genes.TMM.EXPR.matrix - if you use -just_do_TMM_scaling

##################################################################################################################################

cd ${cuwodi}/${outdir}


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 0.05 -C 2 --max_genes_clust 50000 --samples ../samples_described.txt
                                                                                                                                                    
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 1e-2 -C 2 --max_genes_clust 50000 --samples ../samples_described.txt
wait                                                                                                                                                
                                                                                                                                                    
                                                                                                                                                    
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 1e-3 -C 2 --max_genes_clust 50000 --samples ../samples_described.txt
wait                                                                                                                                                
                                                                                                                                                    
                                                                                                                                                    
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 1e-4 -C 2 --max_genes_clust 50000 --samples ../samples_described.txt
wait                                                                                                                                                
                                                                                                                                                    
                                                                                                                                                    
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 1e-5 -C 2 --max_genes_clust 50000 --samples ../samples_described.txt

wait 



###########################################################################################################

# P = 0.05
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 30 -R diffExpr.P0.05_C2.matrix.RData
wait 
mv ${cuwodi}/${outdir}/*_P_30.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C2.matrix.RData.clusters_fixed_P_30

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P0.05_C2.matrix.RData
wait 

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ktree 8 -R diffExpr.P0.05_C2.matrix.RData

mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C2.matrix.RData.clusters_fixed_P_40
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P0.05_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C2.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P0.05_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C2.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P0.05_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C2.matrix.RData.clusters_fixed_P_70
 
wait 


#################################################

# p = 1e-2
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 30 -R diffExpr.P1e-2_C2.matrix.RData
wait 

mv ${cuwodi}/${outdir}/*_P_30.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C2.matrix.RData.clusters_fixed_P_30

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P1e-2_C2.matrix.RData
wait 

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ktree 8 -R diffExpr.P1e-2_C2.matrix.RData


mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C2.matrix.RData.clusters_fixed_P_40
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P1e-2_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C2.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P1e-2_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C2.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P1e-2_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C2.matrix.RData.clusters_fixed_P_70
 
wait 



############################################
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 20 -R diffExpr.P1e-3_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_20.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_20


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 10 -R diffExpr.P1e-3_C2.matrix.RData
wait 



mv ${cuwodi}/${outdir}/*_P_10.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_10


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 30 -R diffExpr.P1e-3_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_30.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_30


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P1e-3_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_40
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P1e-3_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P1e-3_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P1e-3_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_70
wait 


###########################################

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P1e-4_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-4_C2.matrix.RData.clusters_fixed_P_40
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P1e-4_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-4_C2.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P1e-4_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-4_C2.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P1e-4_C2.matrix.RData
wait 


mv *_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-4_C2.matrix.RData.clusters_fixed_P_70
wait 



##########################################################
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P1e-5_C2.matrix.RData

wait 


mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-5_C2.matrix.RData.clusters_fixed_P_40

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P1e-5_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-5_C2.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P1e-5_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-5_C2.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P1e-5_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-5_C2.matrix.RData.clusters_fixed_P_70
wait 


##################################################################

############################################################################################################################################################################################################




outdir=DE_analysis_DESeq2LOG2

cd ${cuwodi}/



/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix *genes.counts.matrix --samples_file samples_described.txt --method DESeq2 --output ${outdir}
wait 


wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --matrix *genes.counts.matrix --lengths genes_feature_lengths.txt


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --lexical_column_ordering --matrix *genes.counts.matrix --samples samples_described.txt --CPM --log2 --min_rowSums 10 --compare_replicates

wait
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --lexical_column_ordering --matrix *genes.counts.matrix --min_rowSums 2 --samples samples_described.txt --log2 --CPM --sample_cor_matrix

wait
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix *genes.counts.matrix -s samples_described.txt --min_rowSums 10 --log2 --CPM --center_rows --prin_comp 3

# does not need gene lenghts for this method
#/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --matrix *genes.counts.matrix --just_do_TMM_scaling


# .genes.counts.matrix.TMM_normalized.FPKM_normalized.FPKM - if you use gene lengths
# *genes.TMM.EXPR.matrix - if you use -just_do_TMM_scaling

##################################################################################################################################

cd ${cuwodi}/${outdir}


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 0.05 -C 2 --max_genes_clust 50000 --samples ../samples_described.txt
                                                                                                                                                    
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 1e-2 -C 2 --max_genes_clust 50000 --samples ../samples_described.txt
wait                                                                                                                                                
                                                                                                                                                    
                                                                                                                                                    
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 1e-3 -C 2 --max_genes_clust 50000 --samples ../samples_described.txt
wait                                                                                                                                                
                                                                                                                                                    
                                                                                                                                                    
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 1e-4 -C 2 --max_genes_clust 50000 --samples ../samples_described.txt
wait                                                                                                                                                
                                                                                                                                                    
                                                                                                                                                    
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 1e-5 -C 2 --max_genes_clust 50000 --samples ../samples_described.txt

wait 



###########################################################################################################

# P = 0.05
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 30 -R diffExpr.P0.05_C2.matrix.RData
wait 
mv ${cuwodi}/${outdir}/*_P_30.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C2.matrix.RData.clusters_fixed_P_30

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P0.05_C2.matrix.RData
wait 

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ktree 8 -R diffExpr.P0.05_C2.matrix.RData

mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C2.matrix.RData.clusters_fixed_P_40
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P0.05_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C2.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P0.05_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C2.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P0.05_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C2.matrix.RData.clusters_fixed_P_70
 
wait 


#################################################

# p = 1e-2
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 30 -R diffExpr.P1e-2_C2.matrix.RData
wait 

mv ${cuwodi}/${outdir}/*_P_30.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C2.matrix.RData.clusters_fixed_P_30

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P1e-2_C2.matrix.RData
wait 

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ktree 8 -R diffExpr.P1e-2_C2.matrix.RData


mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C2.matrix.RData.clusters_fixed_P_40
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P1e-2_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C2.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P1e-2_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C2.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P1e-2_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C2.matrix.RData.clusters_fixed_P_70
 
wait 



############################################
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 20 -R diffExpr.P1e-3_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_20.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_20


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 10 -R diffExpr.P1e-3_C2.matrix.RData
wait 



mv ${cuwodi}/${outdir}/*_P_10.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_10


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 30 -R diffExpr.P1e-3_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_30.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_30


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P1e-3_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_40
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P1e-3_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P1e-3_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P1e-3_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_70
wait 


###########################################

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P1e-4_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-4_C2.matrix.RData.clusters_fixed_P_40
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P1e-4_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-4_C2.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P1e-4_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-4_C2.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P1e-4_C2.matrix.RData
wait 


mv *_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-4_C2.matrix.RData.clusters_fixed_P_70
wait 



##########################################################
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P1e-5_C2.matrix.RData

wait 


mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-5_C2.matrix.RData.clusters_fixed_P_40

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P1e-5_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-5_C2.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P1e-5_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-5_C2.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P1e-5_C2.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-5_C2.matrix.RData.clusters_fixed_P_70
wait 


##################################################################

########################################################################################################################################################################
##############




outdir=DE_analysis_EdgeRLOG1.6

cd ${cuwodi}/


#echo "\tGp_7DPI_rep1\tGp_7DPI_rep2\tGp_14DPI_rep1\tGp_14DPI_rep2\tGp_21DPI_rep1\tGp_21DPI_rep2\tGp_28DPI_rep1\tGp_28DPI_rep2\tGp_35DPI_rep1\tGp_35DPI_rep2\tGp_EGG_rep1\tGp_EGG_rep2\tGp_J2_rep1\tGp_J2_rep2\tGp_J2_rep3\tGp_MALE_rep1\tGp_MALE_rep2" > Gp_genes.counts.matrix
# merge these into a file
#FILES=$(ls -t -v *.counts | tr '\n' ' ')
#awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES | grep -v "__" >> Gp_genes.counts.matrix


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix *genes.counts.matrix --samples_file samples_described.txt --method edgeR --output ${outdir}
wait 


wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --matrix *genes.counts.matrix --lengths genes_feature_lengths.txt


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --lexical_column_ordering --matrix *genes.counts.matrix --samples samples_described.txt --CPM --log2 --min_rowSums 10 --compare_replicates

wait
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --lexical_column_ordering --matrix *genes.counts.matrix --min_rowSums 2 --samples samples_described.txt --log2 --CPM --sample_cor_matrix

wait
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix *genes.counts.matrix -s samples_described.txt --min_rowSums 10 --log2 --CPM --center_rows --prin_comp 3

# does not need gene lenghts for this method
#/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --matrix *genes.counts.matrix --just_do_TMM_scaling


# .genes.counts.matrix.TMM_normalized.FPKM_normalized.FPKM - if you use gene lengths
# *genes.TMM.EXPR.matrix - if you use -just_do_TMM_scaling

##################################################################################################################################

cd ${cuwodi}/${outdir}


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 0.05 -C 1.6 --max_genes_clust 50000 --samples ../samples_described.txt
                                                                                                                                                    
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 1e-2 -C 1.6 --max_genes_clust 50000 --samples ../samples_described.txt
wait                                                                                                                                                
                                                                                                                                                    
                                                                                                                                                    
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 1e-3 -C 1.6 --max_genes_clust 50000 --samples ../samples_described.txt
wait                                                                                                                                                
                                                                                                                                                    
                                                                                                                                                    
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 1e-4 -C 1.6 --max_genes_clust 50000 --samples ../samples_described.txt
wait                                                                                                                                                
                                                                                                                                                    
                                                                                                                                                    
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*genes.TMM.EXPR.matrix -P 1e-5 -C 1.6 --max_genes_clust 50000 --samples ../samples_described.txt

wait 



###########################################################################################################

# P = 0.05
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 30 -R diffExpr.P0.05_C1.6.matrix.RData
wait 
mv ${cuwodi}/${outdir}/*_P_30.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C1.6.matrix.RData.clusters_fixed_P_30

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P0.05_C1.6.matrix.RData
wait 

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ktree 8 -R diffExpr.P0.05_C1.6.matrix.RData

mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C1.6.matrix.RData.clusters_fixed_P_40
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P0.05_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C1.6.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P0.05_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C1.6.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P0.05_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P0.05_C1.6.matrix.RData.clusters_fixed_P_70
 
wait 


#################################################

# p = 1e-2
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 30 -R diffExpr.P1e-2_C1.6.matrix.RData
wait 

mv ${cuwodi}/${outdir}/*_P_30.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C1.6.matrix.RData.clusters_fixed_P_30

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P1e-2_C1.6.matrix.RData
wait 

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ktree 8 -R diffExpr.P1e-2_C1.6.matrix.RData


mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C1.6.matrix.RData.clusters_fixed_P_40
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P1e-2_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C1.6.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P1e-2_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C1.6.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P1e-2_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-2_C1.6.matrix.RData.clusters_fixed_P_70
 
wait 



############################################
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 20 -R diffExpr.P1e-3_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_20.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C1.6.matrix.RData.clusters_fixed_P_20


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 10 -R diffExpr.P1e-3_C1.6.matrix.RData
wait 



mv ${cuwodi}/${outdir}/*_P_10.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C1.6.matrix.RData.clusters_fixed_P_10


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 30 -R diffExpr.P1e-3_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_30.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C1.6.matrix.RData.clusters_fixed_P_30


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P1e-3_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C1.6.matrix.RData.clusters_fixed_P_40
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P1e-3_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C1.6.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P1e-3_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C1.6.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P1e-3_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-3_C1.6.matrix.RData.clusters_fixed_P_70
wait 


###########################################

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P1e-4_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-4_C1.6.matrix.RData.clusters_fixed_P_40
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P1e-4_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-4_C1.6.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P1e-4_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-4_C1.6.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P1e-4_C1.6.matrix.RData
wait 


mv *_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-4_C1.6.matrix.RData.clusters_fixed_P_70
wait 



##########################################################
/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 40 -R diffExpr.P1e-5_C1.6.matrix.RData

wait 


mv ${cuwodi}/${outdir}/*_P_40.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-5_C1.6.matrix.RData.clusters_fixed_P_40

/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R diffExpr.P1e-5_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_50.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-5_C1.6.matrix.RData.clusters_fixed_P_50
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 60 -R diffExpr.P1e-5_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_60.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-5_C1.6.matrix.RData.clusters_fixed_P_60
wait 


/shelf/apps/User_name/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 70 -R diffExpr.P1e-5_C1.6.matrix.RData
wait 


mv ${cuwodi}/${outdir}/*_P_70.heatmap.heatmap.pdf ${cuwodi}/${outdir}/diffExpr.P1e-5_C1.6.matrix.RData.clusters_fixed_P_70
wait 


##################################################################














