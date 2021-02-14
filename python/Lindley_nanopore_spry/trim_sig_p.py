from collections import defaultdict

import os
from sys import stdin,argv




def seq_getter(filename, outfile):
    "this is a function to open up a .xml file blast results, the out put of\
is all the unique hit"
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    sigp = open("all_signal_p.txt", "r")

    gene_sigp = defaultdict()

    for line in sigp:
        if line.startswith("#"): continue
        data = line.split()
        gene = data[0]
        cmax = data[2]
        gene_sigp[gene] = int(cmax) -1
    f= open(outfile, 'w')
    #record = SeqIO.read(filename, "fasta")
    for seq_record in SeqIO.parse(filename, "fasta"):
        #print seq_record.id
        #print 'boomshanka'

        cmax = gene_sigp[seq_record.id]
        if seq_record.id == "GPALN_004734-T1":
            print(len(seq_record.seq), cmax)
        seq_record.seq = seq_record.seq[3*cmax:]
        if seq_record.id == "GPALN_004734-T1":
            print(len(seq_record.seq), cmax)
        
        SeqIO.write(seq_record, f, "fasta")
    f.close()



seq_getter("Newton_SPRYSECS_RBBH_in_LIND_Nanopore.nt.fasta",
           "Newton_SPRYSECS_RBBH_in_LIND_Nanopore.nt.SIGP_removed.fasta")


seq_getter("LIND_nanopore_RBBH_In_newton.nt.fasta",
           "LIND_nanopore_RBBH_In_newton.nt.SIGP_removed.fasta")
