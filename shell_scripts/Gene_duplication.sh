cd /storage/home/users/User_name/newton/gene_models/fun/final_genes/annotations


diamond makedb --in GPAL_newton_v1.0_annot.AA.fasta -d Gp   

diamond blastp -p 4 --sensitive -e 1e-10 -v -q GPAL_newton_v1.0_annot.AA.fasta -d Gp.dmnd -o  Gp_Gp.tab
# individual species Gp
mv Gp_Gp.tab Gp_Gp.blast


 #reformet the GFF files:  - have to modify script if scaffold names dont fit the norm
conda activate python27
python reform_gff.py -m True --gff GPAL_newton_V1.0.gff3 -s Gp -o Gp.gff


# classify gene duplication

/storage/home/users/User_name/shelf_apps/apps/MCScanX/duplicate_gene_classifier  /storage/home/users/User_name/newton/gene_models/fun/final_genes/annotations/Gp/
