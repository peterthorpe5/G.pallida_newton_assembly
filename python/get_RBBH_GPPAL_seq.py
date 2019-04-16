import os
from sys import stdin,argv

def old_to_new(old, old_to_new_dict):
    for line in old:
        if not line.startswith("#"):
            if line.rstrip():
                old = line.split()[1]
                new = line.split()[0]
                new = new.split("-T")[0]
                old_to_new_dict[old] = new
    return old_to_new_dict
        
    

def seq_getter(filename_in, wantedfile, outfile):
    "this is a function to open up a .xml file blast results, the out put of\
is all the unique hit"
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    f= open(outfile, 'w')
    wanted = open(wantedfile, "r")
    old = open("published_vs_newton_RBBH.tab", "r")
    old_to_new_dict = dict()
    old_to_new_dict = old_to_new(old, old_to_new_dict)

    names = wanted.readlines()
    #print names
    wanted_data = [line.replace("\t", "").rstrip("\n") for line in names
              if line.strip() != ""]
    name_set = set([])
    for i in wanted_data:
        if not i.startswith("#"):
            i = i.rstrip()
            if old_to_new_dict.get(i):
                new = old_to_new_dict[i]
                name_set.add(new.rstrip())
            else:
                print("%s does not have a RBBH" % i)
            print(i, new)
    #print wanted_data

    cds_database = SeqIO.index(filename_in, "fasta")
    #record = SeqIO.read(filename, "fasta")
    for i in name_set:
        if "\r\n" in i:
            i = i.replace("\r\n","")
        print i
        if cds_database.get(i):
            seq_record = cds_database[i]
            #print 'boomshanka'
            seq_record.description = ""
            SeqIO.write(seq_record, f, "fasta")
    f.close()
    return True


seq_getter(argv[1],argv[2], argv[3])

#seq_getter('assembly2_scaffolds.fasta',\
           #'scaffold318.fasta')
print 'done'

