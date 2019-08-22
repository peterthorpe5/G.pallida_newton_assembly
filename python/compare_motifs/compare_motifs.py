# title: script to summarise the motifs ound and their closest genes
# authou. Peter Thorpe 2018

from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import sys
from optparse import OptionParser
import datetime
import logging
import logging.handlers
import time


def index_annotation(infasta):
    """function to index the annotaions of the genes of interest"""
    gene_annot = defaultdict(str)
    for seq_record in SeqIO.parse(infasta, "fasta"):
        gene_annot[seq_record.id] = seq_record.description
        gene_annot[seq_record.id.split("-T")[0]] = seq_record.description
    return gene_annot


def parse_summary_file(summary_file, min_motifs=1):
    """function to collect name motif in a specific condition"""
    name_set = set([])
    name_to_annot = dict()
    name_to_motifs = dict()
    gene_upJ2_vs_egg = set([])
    gene_upJ2_vs_14dpi = set([])
    with open(summary_file, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if line.startswith("sample"):
                continue
            if not line.strip():
                continue #  if the last line is blank
            gene, motif_counts, spry, Secreted, \
                annotation, GH, dup, up_J2_v_14dpi, \
                up_J2_v_egg, upregulated = line.split("\t")

            if int(motif_counts) >= min_motifs:
                name_to_motifs[gene.rstrip()] = motif_counts
                if not Secreted == "-":
                    name_set.add(gene.rstrip())
                    
                if not up_J2_v_egg == "-":
                    gene_upJ2_vs_egg.add(gene.rstrip())
                    
                if not up_J2_v_14dpi == "-":
                    gene_upJ2_vs_14dpi.add(gene.rstrip())
                
    return name_set, name_to_annot, name_to_motifs, gene_upJ2_vs_egg, gene_upJ2_vs_14dpi


# run as script
SUBname_set, SUBname_to_annot, SUBname_to_motifs, \
          SUBgene_upJ2_vs_egg, \
          SUBgene_upJ2_vs_14dpi = parse_summary_file("default_upstream_2500_motif1_2500bp_upstream.bed.summary")
DOGname_set, DOGname_to_annot, DOGname_to_motifs, \
          DOGgene_upJ2_vs_egg, \
          DOGgene_upJ2_vs_14dpi = parse_summary_file("DOGBOX_upstream_500_motif_500bp_upstream.bed.summary")

print("Subs secreted, with 1 motif or more: %d" % len(SUBname_set))
print("Subs secreted, with 1 motif or more, UP in J2 versus egg: %d"
      % len(SUBname_set.intersection(SUBgene_upJ2_vs_egg)))
sub_up2 = SUBgene_upJ2_vs_egg.union(SUBgene_upJ2_vs_14dpi)
#print(len(sub_up2), len(SUBgene_upJ2_vs_egg), len(SUBgene_upJ2_vs_14dpi))
#print(len(SUBgene_upJ2_vs_egg))
print("Subs secreted, with 1 motif or more, UP in J2 versus egg and/or up J2 versus 14dpi: %d"
      % len(SUBname_set.intersection(sub_up2)))

print("\n")
###
#dog summary
print("DOGs secreted, with 1 motif or more: %d " % len(DOGname_set))

print("DOG secreted, with 1 motif or more: %d" % len(DOGname_set))
print("DOG secreted, with 1 motif or more, UP in J2 versus egg: %d"
      % len(DOGname_set.intersection(DOGgene_upJ2_vs_egg)))
dog_up2 = DOGgene_upJ2_vs_egg.union(DOGgene_upJ2_vs_14dpi)
#print(len(dog_up2))
#print(len(DOGgene_upJ2_vs_egg))
print("DOG secreted, with 1 motif or more, UP in J2 versus egg and/or up J2 versus 14dpi: %d"
      % len(DOGname_set.intersection(dog_up2)))

print("\n")

# compare dog versus motif summary
print("Sub secreted up J2 vs egg and up J2 versus 14dpi Versus ALL secreted DOGS = %d"
      % (len(sub_up2.intersection(DOGname_set))))


###########################################
### DOGS with 3 motifs or more

DOGname_set, DOGname_to_annot, DOGname_to_motifs, \
          DOGgene_upJ2_vs_egg, \
          DOGgene_upJ2_vs_14dpi = parse_summary_file("DOGBOX_upstream_500_motif_500bp_upstream.bed.summary", 3)

print("Subs secreted, with 1 motif or more: %d" % len(SUBname_set))
print("Subs secreted, with 1 motif or more, UP in J2 versus egg: %d"
      % len(SUBname_set.intersection(SUBgene_upJ2_vs_egg)))
sub_up2 = SUBgene_upJ2_vs_egg.union(SUBgene_upJ2_vs_14dpi)
print(len(sub_up2), len(SUBgene_upJ2_vs_egg), len(SUBgene_upJ2_vs_14dpi))
#print(len(SUBgene_upJ2_vs_egg))
print("Subs secreted, with 1 motif or more, UP in J2 versus egg and/or up J2 versus 14dpi: %d"
      % len(SUBname_set.intersection(sub_up2)))
Final_subvenrtral = SUBname_set.intersection(sub_up2)

print("\n")
###
#dog summary
print("DOGs secreted, with 3 motif or more: %d " % len(DOGname_set))

print("DOG secreted, with 3 motif or more: %d" % len(DOGname_set))
print("DOG secreted, with 3 motif or more, UP in J2 versus egg: %d"
      % len(DOGname_set.intersection(DOGgene_upJ2_vs_egg)))
dog_up2 = DOGgene_upJ2_vs_egg.union(DOGgene_upJ2_vs_14dpi)
#print(len(dog_up2))
#print(len(DOGgene_upJ2_vs_egg))
print("DOG secreted, with 3 motif or more, UP in J2 versus egg and/or up J2 versus 14dpi: %d"
      % len(DOGname_set.intersection(dog_up2)))

print("\n")

# compare dog versus motif summary
print("Sub secreted up J2 vs egg and up J2 versus 14dpi which ARE IN ALL secreted DOGS 3 motifs or more= %d"
      % (len(sub_up2.intersection(DOGname_set))))

print("Sub secreted up J2 vs egg and up J2 versus 14dpi which ARE NOT IN ALL secreted DOGS 3 motifs or more = %d"
      % (len(Final_subvenrtral.difference(DOGname_set))))


print("\n")

annotation = index_annotation("Gpal_newton_annotated.AA.fasta")
Subs_with_DOGS = sub_up2.intersection(DOGname_set)
for i in Subs_with_DOGS:
    #print(i)
    i = i.rstrip()
    annotat = annotation[i.rstrip()]
    print("subvnetral lists with 3 Dogs\t%s\t%s" % (i, annotat))

print("\n")
final_sub_no_dogs = Final_subvenrtral.difference(DOGname_set)
for i in final_sub_no_dogs:
    #print(i)
    i = i.rstrip()
    annotat = annotation[i.rstrip()]
    print("subvnetral lists with NO Dogs\t%s\t%s" % (i, annotat))
    




