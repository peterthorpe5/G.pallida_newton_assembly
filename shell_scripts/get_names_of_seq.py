import os
from sys import stdin,argv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO


def seq_getter(filename1, outfile):
    "this is a function to open up a .xml file blast results, the out put of\
is all the unique hit"

    
    for seq_record in SeqIO.parse(filename1, "fasta"):
        #DNA_seq = seq_record.seq
        #translated_seq = DNA_seq.translate()
        #print translated_seq
        print seq_record.id, "\tlength=\t", len(seq_record.seq)



    f= open(outfile, 'w')
    s_num = 0 
    for seq_record in SeqIO.parse(filename1, "fasta"):
        s_num = s_num + 1
        seq_record.id = "scaffold%d" % s_num
        SeqIO.write(seq_record, f, "fasta")



seq_getter(argv[1],argv[2])

