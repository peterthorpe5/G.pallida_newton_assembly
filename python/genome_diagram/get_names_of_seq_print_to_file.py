
def seq_getter(filename1, outfile):
    "opens up a fasta files and prints all names of sequences to a file."
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Bio import SeqIO
    f_out = open(outfile, "w")
    
    for seq_record in SeqIO.parse(filename1, "fasta"):
        #DNA_seq = seq_record.seq
        #translated_seq = DNA_seq.translate()
        #print translated_seq
        data = seq_record.id + "\t" + str(len(seq_record.seq)) + "\n"
        f_out.write(data)

    return True




seq_getter('Gp_Newton_haplotype1.fasta',\
           'Gp_Newton_haplotype1.seq_len')
