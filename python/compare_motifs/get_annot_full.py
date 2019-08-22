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
    



noneffectors = """GPALN_015211-T1
GPALN_015214-T1
GPALN_010126-T1
GPALN_010127-T1
GPALN_009837-T1
GPALN_002386-T1
GPALN_002387-T1
GPALN_005161-T1
GPALN_010598-T1
GPALN_012099-T1
GPALN_008873-T1
GPALN_014324-T1
GPALN_014327-T1
GPALN_005801-T1
GPALN_004014-T1
GPALN_015654-T1
GPALN_010795-T1
GPALN_013150-T1
GPALN_007748-T1
GPALN_002346-T1
GPALN_004369-T1
GPALN_003369-T1
GPALN_004410-T1
GPALN_004411-T1
GPALN_006782-T1
GPALN_004380-T1
GPALN_006911-T1
GPALN_005160-T1
GPALN_015654-T2
GPALN_010093-T1
GPALN_005162-T1
GPALN_012008-T1
GPALN_010702-T1
GPALN_005042-T1
GPALN_004581-T1
GPALN_001111-T1
GPALN_013350-T1
GPALN_009672-T1
GPALN_015632-T1
GPALN_011858-T1
GPALN_006853-T1
GPALN_009445-T1
GPALN_015107-T1
GPALN_011878-T1
GPALN_010970-T1
GPALN_013383-T1
GPALN_012367-T1
GPALN_012287-T1
GPALN_009458-T1
GPALN_006836-T1
GPALN_010231-T1
GPALN_006824-T1
GPALN_015569-T1
GPALN_006058-T1
GPALN_002293-T1
GPALN_012056-T1
GPALN_013171-T1
GPALN_015278-T1
GPALN_010392-T1
GPALN_008654-T1
GPALN_004893-T1
GPALN_003905-T1
GPALN_003860-T1
GPALN_009443-T1
GPALN_009839-T1
GPALN_002386-T1
GPALN_002387-T1
GPALN_009837-T1
GPALN_003910-T1
GPALN_003908-T1
GPALN_003911-T1
GPALN_006780-T1
GPALN_010789-T1
GPALN_010793-T1
GPALN_013387-T1
GPALN_009668-T1
GPALN_009670-T1
GPALN_012067-T1
GPALN_007796-T1
GPALN_004009-T1
GPALN_006837-T1
GPALN_003787-T1
GPALN_011857-T1
GPALN_004814-T1
GPALN_004010-T1
GPALN_004009-T1
GPALN_007796-T1
GPALN_013387-T1
GPALN_012067-T1
GPALN_009668-T1
GPALN_003787-T1
GPALN_009672-T1
GPALN_015179-T1
GPALN_015182-T1
GPALN_015181-T1
GPALN_003913-T1
GPALN_003946-T1
GPALN_003942-T1
GPALN_003925-T1
GPALN_003943-T1
GPALN_003912-T1
GPALN_015178-T1
GPALN_015188-T1
GPALN_010067-T1
GPALN_003953-T1
GPALN_003949-T1
GPALN_010625-T1
GPALN_003955-T1
GPALN_014539-T1
GPALN_006018-T1
GPALN_016188-T1
GPALN_016185-T1
GPALN_006039-T2
GPALN_016187-T1
GPALN_006039-T1
GPALN_010777-T1
GPALN_014377-T1
GPALN_006037-T1
GPALN_015279-T1
GPALN_007443-T1
GPALN_015297-T1
GPALN_009441-T1
GPALN_009441-T2
GPALN_009441-T3
GPALN_015291-T1
GPALN_002466-T1
GPALN_010067-T1
GPALN_015188-T1
GPALN_003953-T1
GPALN_010625-T1
GPALN_003949-T1
GPALN_003955-T1
GPALN_014539-T1
GPALN_006018-T1
GPALN_002346-T1
GPALN_010795-T1
GPALN_015654-T1
GPALN_004014-T1
GPALN_013150-T1
GPALN_007748-T1
GPALN_004380-T1
GPALN_005160-T1
GPALN_005155-T1
GPALN_003369-T1
GPALN_004369-T1
GPALN_006782-T1
GPALN_002377-T1
GPALN_003891-T1
GPALN_003868-T1
GPALN_003964-T1
GPALN_008468-T1
GPALN_005721-T1
GPALN_004017-T1
GPALN_004018-T1
GPALN_003846-T1
GPALN_014885-T1
GPALN_006778-T1
GPALN_015171-T1
GPALN_002204-T1
GPALN_010793-T1
GPALN_014355-T1
GPALN_014357-T1
GPALN_010414-T1
GPALN_010416-T1
GPALN_010405-T1
GPALN_010412-T1
GPALN_004009-T1
GPALN_007796-T1
GPALN_013387-T1
GPALN_011857-T1
GPALN_009668-T1
GPALN_012067-T1
GPALN_003787-T1
GPALN_009672-T1
GPALN_004010-T1
GPALN_006837-T1
GPALN_009670-T1
GPALN_005903-T1
GPALN_005901-T1
GPALN_015605-T1
GPALN_006742-T1
GPALN_010511-T1
GPALN_010536-T1
GPALN_010603-T1
GPALN_002290-T1
GPALN_014398-T1
GPALN_006828-T1
GPALN_013387-T1
GPALN_009668-T1
GPALN_007796-T1
GPALN_004009-T1
GPALN_003787-T1
GPALN_012067-T1
GPALN_011857-T1
GPALN_009672-T1
GPALN_009670-T1
GPALN_004010-T1
GPALN_006837-T1
GPALN_002593-T1
GPALN_015248-T1
GPALN_003907-T2
GPALN_003907-T1
GPALN_002346-T1
GPALN_010795-T1
GPALN_015654-T1
GPALN_004014-T1
GPALN_007748-T1
GPALN_013150-T1
GPALN_004380-T1
GPALN_005160-T1
GPALN_005155-T1
GPALN_003369-T1
GPALN_004369-T1
GPALN_006782-T1
GPALN_005721-T1
GPALN_003868-T1
GPALN_003964-T1
GPALN_015654-T2
GPALN_003891-T1
GPALN_015179-T1
GPALN_015182-T1
GPALN_015181-T1
GPALN_003913-T1
GPALN_003946-T1
GPALN_003925-T1
GPALN_003942-T1
GPALN_003943-T1
GPALN_003912-T1
GPALN_015178-T1
GPALN_003925-T1
GPALN_003946-T1
GPALN_003943-T1
GPALN_003912-T1
GPALN_003942-T1
GPALN_006035-T1
GPALN_002520-T1
GPALN_010421-T1
GPALN_000143-T1
GPALN_015179-T1
GPALN_015182-T1
GPALN_015181-T1
GPALN_003913-T1
GPALN_003946-T1
GPALN_003942-T1
GPALN_003925-T1
GPALN_003943-T1
GPALN_003912-T1
GPALN_004008-T1
GPALN_003970-T1
GPALN_004493-T1
GPALN_015605-T1
GPALN_006742-T1
GPALN_005901-T1
GPALN_005903-T1
GPALN_009458-T1
GPALN_002300-T1
GPALN_008654-T1
GPALN_011999-T1
GPALN_013342-T1
GPALN_009465-T1
GPALN_010232-T1
GPALN_015107-T1
GPALN_006836-T1
GPALN_015278-T1
GPALN_010361-T1
GPALN_006839-T1
GPALN_015632-T1
GPALN_004981-T1
GPALN_011878-T1
GPALN_004893-T1
GPALN_002293-T1
GPALN_013383-T1
GPALN_015569-T1
GPALN_007441-T1
GPALN_013171-T1
GPALN_011858-T1
GPALN_009445-T1
GPALN_006058-T1
GPALN_009672-T1
GPALN_014865-T1
GPALN_005801-T1
GPALN_014327-T1
GPALN_014324-T1
GPALN_005161-T1
GPALN_010536-T1
GPALN_010511-T1
GPALN_012010-T1
GPALN_002378-T1
GPALN_002322-T1
GPALN_002320-T1
GPALN_004018-T1
GPALN_004017-T1
GPALN_006778-T1
GPALN_008468-T1
GPALN_015171-T1
GPALN_003882-T1
GPALN_014885-T1
GPALN_010171-T1
GPALN_002204-T1
GPALN_003910-T1
GPALN_002377-T1
GPALN_003908-T1
GPALN_003891-T1
GPALN_003846-T1
GPALN_003868-T1
GPALN_007748-T1
GPALN_004369-T1
GPALN_013150-T1
GPALN_003369-T1
GPALN_006782-T1
GPALN_004411-T1
GPALN_015632-T1
GPALN_012056-T1
GPALN_013350-T1
GPALN_013383-T1
GPALN_011878-T1
GPALN_011858-T1
GPALN_009458-T1
GPALN_012287-T1
GPALN_004893-T1
GPALN_013111-T1
GPALN_014271-T1
GPALN_015278-T1
GPALN_010232-T1
GPALN_008654-T1
GPALN_004981-T1
GPALN_015107-T1
GPALN_001382-T1
GPALN_009465-T1
GPALN_015102-T1
GPALN_013171-T1
GPALN_004893-T2
GPALN_010208-T1
GPALN_013342-T1
GPALN_002293-T1
GPALN_010970-T1
GPALN_007648-T1
GPALN_008161-T1
GPALN_010416-T1
GPALN_010414-T1
GPALN_010405-T1
GPALN_010412-T1
GPALN_015654-T1
GPALN_015654-T2
GPALN_010795-T1
GPALN_004014-T1
GPALN_002346-T1
GPALN_003925-T1
GPALN_003946-T1
GPALN_003912-T1
GPALN_003943-T1
GPALN_003942-T1
GPALN_003913-T1
GPALN_015182-T1
GPALN_015179-T1
GPALN_015181-T1
GPALN_010171-T1
GPALN_015171-T1
GPALN_003882-T1
GPALN_014885-T1
GPALN_008468-T1
GPALN_010795-T1
GPALN_003891-T1
GPALN_012010-T1
GPALN_004018-T1
GPALN_015654-T1
GPALN_002204-T1
GPALN_004017-T1
GPALN_002191-T1
GPALN_006778-T1
GPALN_003868-T1
GPALN_002378-T1
GPALN_003964-T1
GPALN_013150-T1
GPALN_002377-T1
GPALN_003846-T1
GPALN_002322-T1
GPALN_002202-T1
GPALN_004380-T1
GPALN_003908-T1
GPALN_005160-T1
GPALN_002386-T1
GPALN_002387-T1
GPALN_009837-T1
GPALN_005161-T1
GPALN_010598-T1
GPALN_012099-T1
GPALN_010067-T1
GPALN_015188-T1
GPALN_003953-T1
GPALN_003949-T1
GPALN_010625-T1
GPALN_014539-T1
GPALN_003955-T1
GPALN_006018-T1
GPALN_007647-T1
GPALN_015425-T1
GPALN_014034-T1
GPALN_016298-T1
GPALN_005901-T1
GPALN_005903-T1
GPALN_015605-T1
GPALN_006742-T1
GPALN_004712-T1
GPALN_011865-T1
GPALN_011865-T2
GPALN_012594-T1
GPALN_012594-T2
GPALN_013144-T1
GPALN_004470-T1
GPALN_010536-T1
GPALN_010603-T1
GPALN_010511-T1
GPALN_002377-T1
GPALN_002379-T1
GPALN_010555-T1
GPALN_003871-T1
GPALN_014713-T1
GPALN_011865-T1
GPALN_011865-T2
GPALN_014576-T1
GPALN_015605-T1
GPALN_006742-T1
GPALN_005901-T1
GPALN_005903-T1
GPALN_004018-T1
GPALN_004017-T1
GPALN_006778-T1
GPALN_014885-T1
GPALN_002378-T1
GPALN_016343-T1
""".split()

print("GGEEETTTIIINGG ANNOT\n\n\n")
for i in noneffectors:
    print(annotation[i.rstrip()])


