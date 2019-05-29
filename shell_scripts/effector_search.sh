

GRSO effectors from here: https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-0985-1/MediaObjects/13059_2016_985_MOESM12_ESM.xlsx
# strict search 
(python27) pjt6@phylo:~/newton/final_genes/funaannot/update_results > diamond blastp --query GROS_High_con_effectors.fasta -d Gp_new.dmnd -e 1e-80 --outfmt 6 --out \
 Gp_newt_vs_GRSO_HC_effectors.1e80.tab

# loose search 
diamond blastp --query GROS_High_con_effectors.fasta -d Gp_new.dmnd -e 1e-20 --outfmt 6 --out Gp_newt_vs_GRSO_HC_effectors.1e20.tab

cut -f2 Gp_newt_vs_GRSO_HC_effectors.1e80.tab > Gp_newt_HC_GROS_1e80.hits


cut -f2 Gp_newt_vs_GRSO_HC_effectors.1e20.tab > Gp_newt_HC_GROS_1e20.hits
