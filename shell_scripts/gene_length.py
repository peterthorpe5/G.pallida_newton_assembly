#title: script to get names and length of gene for trinity DE analysis

""" Why?
The genes and lengths are required for normalisation.

the output should be like this:

gene_id	length	effective_length
Rpadi_10000|c0_g1	1826.75	1635.59
Rpadi_10000|c1_g1	296.00	106.94
Rpadi_10000|c3_g1	426.00	235.01
Rpadi_10000|c4_g1	303.00	113.66
Rpadi_10003|c0_g1	2946.69	2755.51

For the predicted genes in the genome, the length == effective length.
Therefore these are
ouputted as the same value
"""
#########################################################################
#imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser

####################################################################################

def seq_getter(filename_in, outfile):
    "opens up a fasta files and prints all names of sequences to a file."
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    f= open(outfile, 'w')
    f.write("gene_id	length	effective_length\n")
    
    for seq_record in SeqIO.parse(filename_in, "fasta"):
        data_formatted = "%s\t%d\t%d\n" %(seq_record.id,
                                    len(seq_record.seq), len(seq_record))
        f.write(data_formatted)

    return True

###############################################################################


if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ python gene_length.py -i gene.fasta -o gene.lengths (default = genes_feature_lengths.txt)

script to gene the gen names and lengths for
trinity normalisation
"""

parser = OptionParser(usage=usage)

parser.add_option("-o", "--out", dest="outfile", default="genes_feature_lengths.txt",
                  help="outfile_name")

parser.add_option("-i", "--in", dest="filename_in", default=None,
                  help="in_file.fasta",
                  metavar="FILE")






(options, args) = parser.parse_args()

outfile = options.outfile
filename_in = options.filename_in

# call the function
seq_getter(filename_in, outfile)

print 'done'
