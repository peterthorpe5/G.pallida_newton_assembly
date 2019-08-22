cd /storage/home/users/pjt6/synteny
cat hg_hg.fa li_li.fa > hg_li.fa
diamond makedb --in hg_li.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q hg_li.fa -d temp.dmnd --outfmt 6 -o ./hg_li/hg_li.blast --max-target-seqs  5 
rm temp*
cat hg_hg.fa gr_gr.fa > hg_gr.fa
diamond makedb --in hg_gr.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q hg_gr.fa -d temp.dmnd --outfmt 6 -o ./hg_gr/hg_gr.blast --max-target-seqs  5 
rm temp*
cat hg_hg.fa mi_mi.fa > hg_mi.fa
diamond makedb --in hg_mi.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q hg_mi.fa -d temp.dmnd --outfmt 6 -o ./hg_mi/hg_mi.blast --max-target-seqs  5 
rm temp*
cat hg_hg.fa gp_gp.fa > hg_gp.fa
diamond makedb --in hg_gp.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q hg_gp.fa -d temp.dmnd --outfmt 6 -o ./hg_gp/hg_gp.blast --max-target-seqs  5 
rm temp*
cat li_li.fa hg_hg.fa > li_hg.fa
diamond makedb --in li_hg.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q li_hg.fa -d temp.dmnd --outfmt 6 -o ./li_hg/li_hg.blast --max-target-seqs  5 
rm temp*
cat li_li.fa gr_gr.fa > li_gr.fa
diamond makedb --in li_gr.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q li_gr.fa -d temp.dmnd --outfmt 6 -o ./li_gr/li_gr.blast --max-target-seqs  5 
rm temp*
cat li_li.fa mi_mi.fa > li_mi.fa
diamond makedb --in li_mi.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q li_mi.fa -d temp.dmnd --outfmt 6 -o ./li_mi/li_mi.blast --max-target-seqs  5 
rm temp*
cat li_li.fa gp_gp.fa > li_gp.fa
diamond makedb --in li_gp.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q li_gp.fa -d temp.dmnd --outfmt 6 -o ./li_gp/li_gp.blast --max-target-seqs  5 
rm temp*
cat gr_gr.fa hg_hg.fa > gr_hg.fa
diamond makedb --in gr_hg.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q gr_hg.fa -d temp.dmnd --outfmt 6 -o ./gr_hg/gr_hg.blast --max-target-seqs  5 
rm temp*
cat gr_gr.fa li_li.fa > gr_li.fa
diamond makedb --in gr_li.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q gr_li.fa -d temp.dmnd --outfmt 6 -o ./gr_li/gr_li.blast --max-target-seqs  5 
rm temp*
cat gr_gr.fa mi_mi.fa > gr_mi.fa
diamond makedb --in gr_mi.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q gr_mi.fa -d temp.dmnd --outfmt 6 -o ./gr_mi/gr_mi.blast --max-target-seqs  5 
rm temp*
cat gr_gr.fa gp_gp.fa > gr_gp.fa
diamond makedb --in gr_gp.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q gr_gp.fa -d temp.dmnd --outfmt 6 -o ./gr_gp/gr_gp.blast --max-target-seqs  5 
rm temp*
cat mi_mi.fa hg_hg.fa > mi_hg.fa
diamond makedb --in mi_hg.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q mi_hg.fa -d temp.dmnd --outfmt 6 -o ./mi_hg/mi_hg.blast --max-target-seqs  5 
rm temp*
cat mi_mi.fa li_li.fa > mi_li.fa
diamond makedb --in mi_li.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q mi_li.fa -d temp.dmnd --outfmt 6 -o ./mi_li/mi_li.blast --max-target-seqs  5 
rm temp*
cat mi_mi.fa gr_gr.fa > mi_gr.fa
diamond makedb --in mi_gr.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q mi_gr.fa -d temp.dmnd --outfmt 6 -o ./mi_gr/mi_gr.blast --max-target-seqs  5 
rm temp*
cat mi_mi.fa gp_gp.fa > mi_gp.fa
diamond makedb --in mi_gp.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q mi_gp.fa -d temp.dmnd --outfmt 6 -o ./mi_gp/mi_gp.blast --max-target-seqs  5 
rm temp*
cat gp_gp.fa hg_hg.fa > gp_hg.fa
diamond makedb --in gp_hg.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q gp_hg.fa -d temp.dmnd --outfmt 6 -o ./gp_hg/gp_hg.blast --max-target-seqs  5 
rm temp*
cat gp_gp.fa li_li.fa > gp_li.fa
diamond makedb --in gp_li.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q gp_li.fa -d temp.dmnd --outfmt 6 -o ./gp_li/gp_li.blast --max-target-seqs  5 
rm temp*
cat gp_gp.fa gr_gr.fa > gp_gr.fa
diamond makedb --in gp_gr.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q gp_gr.fa -d temp.dmnd --outfmt 6 -o ./gp_gr/gp_gr.blast --max-target-seqs  5 
rm temp*
cat gp_gp.fa mi_mi.fa > gp_mi.fa
diamond makedb --in gp_mi.fa -d temp 
diamond blastp -p 12 --more-sensitive -e 1e-5 -v -q gp_mi.fa -d temp.dmnd --outfmt 6 -o ./gp_mi/gp_mi.blast --max-target-seqs  5 
rm temp*

##############
#mkdir gp_gr 
#cat gp_gp.gff gr_gr.gff > gp_gr.gff
#cat gp_gp.gff gr_gr.gff > ./gp_gr/gp_gr.gff
#mkdir gp_hg 
#cat gp_gp.gff hg_hg.gff > gp_hg.gff
#cat gp_gp.gff hg_hg.gff > ./gp_hg/gp_hg.gff
#mkdir gp_li 
#cat gp_gp.gff li_li.gff > gp_li.gff
#cat gp_gp.gff li_li.gff > ./gp_li/gp_li.gff
#mkdir gp_mi 
#cat gp_gp.gff mi_mi.gff > gp_mi.gff
#cat gp_gp.gff mi_mi.gff > ./gp_mi/gp_mi.gff
#mkdir gr_gp 
#cat gr_gr.gff gp_gp.gff > gr_gp.gff
#cat gr_gr.gff gp_gp.gff > ./gr_gp/gr_gp.gff
#mkdir gr_hg 
#cat gr_gr.gff hg_hg.gff > gr_hg.gff
#cat gr_gr.gff hg_hg.gff > ./gr_hg/gr_hg.gff
#mkdir gr_li 
#cat gr_gr.gff li_li.gff > gr_li.gff
#cat gr_gr.gff li_li.gff > ./gr_li/gr_li.gff
#mkdir gr_mi 
#cat gr_gr.gff mi_mi.gff > gr_mi.gff
#cat gr_gr.gff mi_mi.gff > ./gr_mi/gr_mi.gff
#mkdir hg_gp 
#cat hg_hg.gff gp_gp.gff > hg_gp.gff
#cat hg_hg.gff gp_gp.gff > ./hg_gp/hg_gp.gff
#mkdir hg_gr 
#cat hg_hg.gff gr_gr.gff > hg_gr.gff
#cat hg_hg.gff gr_gr.gff > ./hg_gr/hg_gr.gff
#mkdir hg_li 
#cat hg_hg.gff li_li.gff > hg_li.gff
#cat hg_hg.gff li_li.gff > ./hg_li/hg_li.gff
#mkdir hg_mi 
#cat hg_hg.gff mi_mi.gff > hg_mi.gff
#cat hg_hg.gff mi_mi.gff > ./hg_mi/hg_mi.gff
#mkdir li_gp 
#cat li_li.gff gp_gp.gff > li_gp.gff
#cat li_li.gff gp_gp.gff > ./li_gp/li_gp.gff
#mkdir li_gr 
#cat li_li.gff gr_gr.gff > li_gr.gff
#cat li_li.gff gr_gr.gff > ./li_gr/li_gr.gff
#mkdir li_hg 
#cat li_li.gff hg_hg.gff > li_hg.gff
#cat li_li.gff hg_hg.gff > ./li_hg/li_hg.gff
#mkdir li_mi 
#cat li_li.gff mi_mi.gff > li_mi.gff
#cat li_li.gff mi_mi.gff > ./li_mi/li_mi.gff
#mkdir mi_gp 
#cat mi_mi.gff gp_gp.gff > mi_gp.gff
#cat mi_mi.gff gp_gp.gff > ./mi_gp/mi_gp.gff
#mkdir mi_gr 
#cat mi_mi.gff gr_gr.gff > mi_gr.gff
#cat mi_mi.gff gr_gr.gff > ./mi_gr/mi_gr.gff
#mkdir mi_hg 
#cat mi_mi.gff hg_hg.gff > mi_hg.gff
#cat mi_mi.gff hg_hg.gff > ./mi_hg/mi_hg.gff
#mkdir mi_li 
#cat mi_mi.gff li_li.gff > mi_li.gff
#cat mi_mi.gff li_li.gff > ./mi_li/mi_li.gff