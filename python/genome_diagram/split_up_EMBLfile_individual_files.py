#title: this script splits up the huge EMBL file made by the perl scripts
 # into one file per contig

#why? so these can be drawn with genomediagram


import os
from sys import stdin,argv

# make a folder for the out file. There could be thousands

file_name = 'test.txt'
script_dir = os.path.dirname(os.path.abspath(__file__))
dest_dir = os.path.join(script_dir, 'individual_files')
try:
    os.makedirs(dest_dir)
except OSError:
    print("already exists")


def scaff_to_len(in_file):
    """open file sacff\tlen
    save to dictionary"""
    f_in = open(in_file, "r")
    data = f_in.readlines()
    scaff_to_len_dict = dict()
    for line in data:
        scaff, length = line.split("\t")
        scaff_to_len_dict[scaff.rstrip()] = length.rstrip()
    return scaff_to_len_dict 

    
def emble_file_splitter(embl_file, scaff_to_len_dict):
    """input is the huge embl file.
    out are individual files fo each contig"""
    f_in = open(embl_file, "r")
    data = f_in.readlines()
    f_out = False
    for line in data:
        if line.startswith("ID   "):
            new_file_name = line.split("ID   ")[1]
            new_file_name = new_file_name.split("|")[0]
            out_name = "./individual_files/%s.embl" % (new_file_name.rstrip("\n"))
            if f_out:
                f_out.close()
            # this name format is required for the genome diagram part.
            # need len of scaff at end tooo
            length = scaff_to_len_dict[new_file_name.rstrip()]
            title = line.rstrip("\n") + "; SV ; ; ; ; ; %s BP.\n" % (str(length))
            f_out = open(out_name, "w")
            f_out.write(title)
        else:
            f_out.write(line)


       
print("usage: python plit_up_EMBL...py file.EMBL")

scaff_to_len_dict = scaff_to_len('Gp_Newton_haplotype1.seq_len')
emble_file_splitter("Gpal_newton.embl", scaff_to_len_dict)
#emble_file_splitter(argv[1])
