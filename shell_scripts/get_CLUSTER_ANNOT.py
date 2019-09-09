import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

handle_name = "Cluster_catorgory_summary.out" 
handle = open(handle_name,"w")
count = 0

cds_database = SeqIO.index("all.annotated.fasta", "fasta")
f_in = open("Orthogroups.txt", "r")
for line in f_in:
    count = count + 1
    data = line.split()
    elements = data[0:]
    annot = set([])
    for gene in elements:
        gene = gene.rstrip()
        if gene in cds_database:
            seq_record = cds_database[gene]
            annot.add(seq_record.description)
    out_anno = ""
    out_count = 0 
    for i in annot:
        #print(i)
        out_count = out_count + 1
        if out_count < 4:
            #print(out_anno)
            out_anno = out_anno + " " + i.rstrip()
    out_fmt = "%d\t%s\n" % (count, out_anno)
    print(out_fmt)
    handle.write(out_fmt)

        
handle.close()
