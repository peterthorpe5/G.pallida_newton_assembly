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
    
    
def parse_DE_file(DE_file):
    """function to collect name DE in a specific condition"""
    DE_name_set = set([])
    with open(DE_file, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if line.startswith("sample"):
                continue
            if not line.strip():
                continue #  if the last line is blank
            gene = line.split()[0]
            DE_name_set.add(gene.rstrip())
    return DE_name_set


def parse_file(infile):
    """function to count the number of times a gene is found frome
    scnatofis. And puill out the annotation"""
    gene_motif_counter = defaultdict(int)
    with open(infile, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if not line.strip():
                continue #  if the last line is blank
            gene = line.rstrip()
            gene_motif_counter[gene] += 1
    return gene_motif_counter
    

def parse_GH_file(cayzymes):
    """function to collect the GH domain info"""
    gene_GH_domian = defaultdict(str)
    with open(cayzymes, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if not line.strip():
                continue #  if the last line is blank
            gene = line.split()[0]
            GH = line.split()[2]
            outinfo = GH + "; "
            gene_GH_domian[gene.split("-T")[0]] += outinfo
    return gene_GH_domian


def parse_duplication_file(dup):
    """function to collect the GH domain info"""
    gene_duplication = defaultdict(str)
    with open(dup, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if not line.strip():
                continue #  if the last line is blank
            gene = line.split()[0]
            dup = line.split()[1]
            outinfo = dup
            gene_duplication[gene.split("-T")[0]] += outinfo
    return gene_duplication
 
def add_upreg_info(upregulated, condition):
    if upregulated == "-":
        upregulated = condition
        return upregulated
    else:
        upregulated = upregulated + "; " + condition
        return upregulated

def print_out_data(gene_annot, gene_motif_counter, 
                  SPRYSECS, sig_and_phobius, 
                  gene_GH_domain, gene_duplication,
                  dpi14_vs_Gp_dpi21_Gp_dpi14_UP, dpi14_vs_Gp_dpi21_Gp_dpi21_UP, 
                  dpi14_vs_Gp_dpi28_Gp_dpi14_UP, dpi14_vs_Gp_dpi28_Gp_dpi28_UP, 
                  dpi14_vs_Gp_dpi35_Gp_dpi14_UP, dpi14_vs_Gp_dpi35_Gp_dpi35_UP, 
                  dpi14_vs_Gp_dpi7_Gp_dpi14_UP, dpi14_vs_Gp_dpi7_Gp_dpi7_UP, 
                  dpi14_vs_Gp_EGG_Gp_dpi14_UP, dpi14_vs_Gp_EGG_Gp_EGG_UP, 
                  dpi14_vs_Gp_J2_Gp_dpi14_UP, dpi14_vs_Gp_J2_Gp_J2_UP, 
                  dpi14_vs_Gp_MALE_Gp_dpi14_UP, dpi14_vs_Gp_MALE_Gp_MALE_UP, 
                  dpi21_vs_Gp_dpi28_Gp_dpi21_UP, dpi21_vs_Gp_dpi28_Gp_dpi28_UP, 
                  dpi21_vs_Gp_dpi35_Gp_dpi21_UP, dpi21_vs_Gp_dpi35_Gp_dpi35_UP,
                  dpi21_vs_Gp_dpi7_Gp_dpi21_UP, dpi21_vs_Gp_dpi7_Gp_dpi7_UP, 
                  dpi21_vs_Gp_EGG_Gp_dpi21_UP, dpi21_vs_Gp_EGG_Gp_EGG_UP, 
                  dpi21_vs_Gp_J2_Gp_dpi21_UP, dpi21_vs_Gp_J2_Gp_J2_UP, 
                  dpi21_vs_Gp_MALE_Gp_dpi21_UP, dpi21_vs_Gp_MALE_Gp_MALE_UP, 
                  dpi28_vs_Gp_dpi35_Gp_dpi28_UP, dpi28_vs_Gp_dpi35_Gp_dpi35_UP, 
                  dpi28_vs_Gp_dpi7_Gp_dpi28_UP, dpi28_vs_Gp_dpi7_Gp_dpi7_UP, 
                  dpi28_vs_Gp_EGG_Gp_dpi28_UP, dpi28_vs_Gp_EGG_Gp_EGG_UP, 
                  dpi28_vs_Gp_J2_Gp_dpi28_UP, dpi28_vs_Gp_J2_Gp_J2_UP, 
                  dpi28_vs_Gp_MALE_Gp_dpi28_UP, dpi28_vs_Gp_MALE_Gp_MALE_UP, 
                  dpi35_vs_Gp_dpi7_Gp_dpi35_UP, dpi35_vs_Gp_dpi7_Gp_dpi7_UP, 
                  dpi35_vs_Gp_EGG_Gp_dpi35_UP, dpi35_vs_Gp_EGG_Gp_EGG_UP, 
                  dpi35_vs_Gp_J2_Gp_dpi35_UP, dpi35_vs_Gp_J2_Gp_J2_UP, 
                  dpi35_vs_Gp_MALE_Gp_dpi35_UP, dpi35_vs_Gp_MALE_Gp_MALE_UP, 
                  dpi7_vs_Gp_EGG_Gp_dpi7_UP, dpi7_vs_Gp_EGG_Gp_EGG_UP, 
                  dpi7_vs_Gp_J2_Gp_dpi7_UP, dpi7_vs_Gp_J2_Gp_J2_UP, 
                  dpi7_vs_Gp_MALE_Gp_dpi7_UP, dpi7_vs_Gp_MALE_Gp_MALE_UP, 
                  EGG_vs_Gp_J2_Gp_EGG_UP, EGG_vs_Gp_J2_Gp_J2_UP, 
                  EGG_vs_Gp_MALE_Gp_EGG_UP, EGG_vs_Gp_MALE_Gp_MALE_UP, 
                  J2_vs_Gp_MALE_Gp_J2_UP, J2_vs_Gp_MALE_Gp_MALE_UP,
                  out):
    """funct to print out the number of times a motif occurs
    and the gene annot"""
    f_out = open(out, "w")
    f_out.write("#gene\tnumber_of_DOGBOX_motif_upstream\tnumber_of_SUBBOX_motif_upstream\tSPRYSEC(Y/N)\tSecreted(Y/n)\tannotation\tGH_domain\tduplication_type\tup_J2_v_14dpi\tup_J2_v_egg\tUp_14_vs_J2\t Up_14_vs_EGG\tup_male_vs_dpi14\tupregulated_LOG2_FDR:p<0.001\n")
    gene_to_dogboxes = parse_duplication_file("./data/secreted_DOG_BOX.txt")
    gene_to_subboxes = parse_duplication_file("./data/secreted_SUBVENTRAL_box.txt")
    for gene, motif_counts in sorted(gene_motif_counter.items()):
        spry = "-"
        Secreted = "-"
        upregulated = "-"
        up_J2_v_14dpi = "-"
        up_J2_v_egg = "-"
        Up_dpi14_vs_Gp_J2 = "-"
        Up_dpi14_vs_Gp_EGG = "-"
        up_male_vs_dpi14 = "-"
        dogs = "-"
        subs = "-"
        annotation = gene_annot[gene]
        gene_temp = gene + "-T1"
        if gene in gene_to_dogboxes:
            dogs = gene_to_dogboxes[gene]
        if gene in gene_to_subboxes:
            subs = gene_to_subboxes[gene]
        if gene_temp in SPRYSECS:
            spry = "YES"
        if gene_temp in sig_and_phobius:
            Secreted = "YES"
        GH = gene_GH_domain[gene.rstrip()]
        if GH == "":
            GH = "-"
        if gene in dpi14_vs_Gp_dpi21_Gp_dpi14_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_dpi21_Gp_dpi14_UP')

        if gene in dpi14_vs_Gp_dpi21_Gp_dpi21_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_dpi21_Gp_dpi21_UP')

        if gene in dpi14_vs_Gp_dpi28_Gp_dpi14_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_dpi28_Gp_dpi14_UP')

        if gene in dpi14_vs_Gp_dpi28_Gp_dpi28_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_dpi28_Gp_dpi28_UP')

        if gene in dpi14_vs_Gp_dpi35_Gp_dpi14_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_dpi35_Gp_dpi14_UP')

        if gene in dpi14_vs_Gp_dpi35_Gp_dpi35_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_dpi35_Gp_dpi35_UP')

        if gene in dpi14_vs_Gp_dpi7_Gp_dpi14_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_dpi7_Gp_dpi14_UP')

        if gene in dpi14_vs_Gp_dpi7_Gp_dpi7_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_dpi7_Gp_dpi7_UP')

        if gene in dpi14_vs_Gp_EGG_Gp_dpi14_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_EGG_Gp_dpi14_UP')
            Up_dpi14_vs_Gp_EGG = "Up_dpi14_vs_EGG"

        if gene in dpi14_vs_Gp_EGG_Gp_EGG_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_EGG_Gp_EGG_UP')

        if gene in dpi14_vs_Gp_J2_Gp_dpi14_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_J2_Gp_dpi14_UP')
            Up_dpi14_vs_Gp_J2 = "Up_14dpi_vs_J2"

        if gene in dpi14_vs_Gp_J2_Gp_J2_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_J2_Gp_J2_UP')
            up_J2_v_14dpi = "up_J2_vs_14dpi"

        if gene in dpi14_vs_Gp_MALE_Gp_dpi14_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_MALE_Gp_dpi14_UP')
            up_male_vs_dpi14 = "up_male_vs_dpi14"

        if gene in dpi14_vs_Gp_MALE_Gp_MALE_UP:
            upregulated = add_upreg_info(upregulated, 'dpi14_vs_Gp_MALE_Gp_MALE_UP')

        if gene in dpi21_vs_Gp_dpi28_Gp_dpi21_UP:
            upregulated = add_upreg_info(upregulated, 'dpi21_vs_Gp_dpi28_Gp_dpi21_UP')

        if gene in dpi21_vs_Gp_dpi28_Gp_dpi28_UP:
            upregulated = add_upreg_info(upregulated, 'dpi21_vs_Gp_dpi28_Gp_dpi28_UP')

        if gene in dpi21_vs_Gp_dpi35_Gp_dpi21_UP:
            upregulated = add_upreg_info(upregulated, 'dpi21_vs_Gp_dpi35_Gp_dpi21_UP')

        if gene in dpi21_vs_Gp_dpi35_Gp_dpi35_UP:
            upregulated = add_upreg_info(upregulated, 'dpi21_vs_Gp_dpi35_Gp_dpi35_UP')

        if gene in dpi21_vs_Gp_dpi7_Gp_dpi21_UP:
            upregulated = add_upreg_info(upregulated, 'dpi21_vs_Gp_dpi7_Gp_dpi21_UP')

        if gene in dpi21_vs_Gp_dpi7_Gp_dpi7_UP:
            upregulated = add_upreg_info(upregulated, 'dpi21_vs_Gp_dpi7_Gp_dpi7_UP')

        if gene in dpi21_vs_Gp_EGG_Gp_dpi21_UP:
            upregulated = add_upreg_info(upregulated, 'dpi21_vs_Gp_EGG_Gp_dpi21_UP')

        if gene in dpi21_vs_Gp_EGG_Gp_EGG_UP:
            upregulated = add_upreg_info(upregulated, 'dpi21_vs_Gp_EGG_Gp_EGG_UP')

        if gene in dpi21_vs_Gp_J2_Gp_dpi21_UP:
            upregulated = add_upreg_info(upregulated, 'dpi21_vs_Gp_J2_Gp_dpi21_UP')

        if gene in dpi21_vs_Gp_J2_Gp_J2_UP:
            upregulated = add_upreg_info(upregulated, 'dpi21_vs_Gp_J2_Gp_J2_UP')

        if gene in dpi21_vs_Gp_MALE_Gp_dpi21_UP:
            upregulated = add_upreg_info(upregulated, 'dpi21_vs_Gp_MALE_Gp_dpi21_UP')

        if gene in dpi21_vs_Gp_MALE_Gp_MALE_UP:
            upregulated = add_upreg_info(upregulated, 'dpi21_vs_Gp_MALE_Gp_MALE_UP')

        if gene in dpi28_vs_Gp_dpi35_Gp_dpi28_UP:
            upregulated = add_upreg_info(upregulated, 'dpi28_vs_Gp_dpi35_Gp_dpi28_UP')

        if gene in dpi28_vs_Gp_dpi35_Gp_dpi35_UP:
            upregulated = add_upreg_info(upregulated, 'dpi28_vs_Gp_dpi35_Gp_dpi35_UP')

        if gene in dpi28_vs_Gp_dpi7_Gp_dpi28_UP:
            upregulated = add_upreg_info(upregulated, 'dpi28_vs_Gp_dpi7_Gp_dpi28_UP')

        if gene in dpi28_vs_Gp_dpi7_Gp_dpi7_UP:
            upregulated = add_upreg_info(upregulated, 'dpi28_vs_Gp_dpi7_Gp_dpi7_UP')

        if gene in dpi28_vs_Gp_EGG_Gp_dpi28_UP:
            upregulated = add_upreg_info(upregulated, 'dpi28_vs_Gp_EGG_Gp_dpi28_UP')

        if gene in dpi28_vs_Gp_EGG_Gp_EGG_UP:
            upregulated = add_upreg_info(upregulated, 'dpi28_vs_Gp_EGG_Gp_EGG_UP')

        if gene in dpi28_vs_Gp_J2_Gp_dpi28_UP:
            upregulated = add_upreg_info(upregulated, 'dpi28_vs_Gp_J2_Gp_dpi28_UP')

        if gene in dpi28_vs_Gp_J2_Gp_J2_UP:
            upregulated = add_upreg_info(upregulated, 'dpi28_vs_Gp_J2_Gp_J2_UP')

        if gene in dpi28_vs_Gp_MALE_Gp_dpi28_UP:
            upregulated = add_upreg_info(upregulated, 'dpi28_vs_Gp_MALE_Gp_dpi28_UP')

        if gene in dpi28_vs_Gp_MALE_Gp_MALE_UP:
            upregulated = add_upreg_info(upregulated, 'dpi28_vs_Gp_MALE_Gp_MALE_UP')

        if gene in dpi35_vs_Gp_dpi7_Gp_dpi35_UP:
            upregulated = add_upreg_info(upregulated, 'dpi35_vs_Gp_dpi7_Gp_dpi35_UP')

        if gene in dpi35_vs_Gp_dpi7_Gp_dpi7_UP:
            upregulated = add_upreg_info(upregulated, 'dpi35_vs_Gp_dpi7_Gp_dpi7_UP')

        if gene in dpi35_vs_Gp_EGG_Gp_dpi35_UP:
            upregulated = add_upreg_info(upregulated, 'dpi35_vs_Gp_EGG_Gp_dpi35_UP')

        if gene in dpi35_vs_Gp_EGG_Gp_EGG_UP:
            upregulated = add_upreg_info(upregulated, 'dpi35_vs_Gp_EGG_Gp_EGG_UP')

        if gene in dpi35_vs_Gp_J2_Gp_dpi35_UP:
            upregulated = add_upreg_info(upregulated, 'dpi35_vs_Gp_J2_Gp_dpi35_UP')

        if gene in dpi35_vs_Gp_J2_Gp_J2_UP:
            upregulated = add_upreg_info(upregulated, 'dpi35_vs_Gp_J2_Gp_J2_UP')

        if gene in dpi35_vs_Gp_MALE_Gp_dpi35_UP:
            upregulated = add_upreg_info(upregulated, 'dpi35_vs_Gp_MALE_Gp_dpi35_UP')

        if gene in dpi35_vs_Gp_MALE_Gp_MALE_UP:
            upregulated = add_upreg_info(upregulated, 'dpi35_vs_Gp_MALE_Gp_MALE_UP')

        if gene in dpi7_vs_Gp_EGG_Gp_dpi7_UP:
            upregulated = add_upreg_info(upregulated, 'dpi7_vs_Gp_EGG_Gp_dpi7_UP')

        if gene in dpi7_vs_Gp_EGG_Gp_EGG_UP:
            upregulated = add_upreg_info(upregulated, 'dpi7_vs_Gp_EGG_Gp_EGG_UP')

        if gene in dpi7_vs_Gp_J2_Gp_dpi7_UP:
            upregulated = add_upreg_info(upregulated, 'dpi7_vs_Gp_J2_Gp_dpi7_UP')

        if gene in dpi7_vs_Gp_J2_Gp_J2_UP:
            upregulated = add_upreg_info(upregulated, 'dpi7_vs_Gp_J2_Gp_J2_UP')

        if gene in dpi7_vs_Gp_MALE_Gp_dpi7_UP:
            upregulated = add_upreg_info(upregulated, 'dpi7_vs_Gp_MALE_Gp_dpi7_UP')

        if gene in dpi7_vs_Gp_MALE_Gp_MALE_UP:
            upregulated = add_upreg_info(upregulated, 'dpi7_vs_Gp_MALE_Gp_MALE_UP')

        if gene in EGG_vs_Gp_J2_Gp_EGG_UP:
            upregulated = add_upreg_info(upregulated, 'EGG_vs_Gp_J2_Gp_EGG_UP')

        if gene in EGG_vs_Gp_J2_Gp_J2_UP:
            upregulated = add_upreg_info(upregulated, 'EGG_vs_Gp_J2_Gp_J2_UP')
            up_J2_v_egg = "up_j2_v_egg"

        if gene in EGG_vs_Gp_MALE_Gp_EGG_UP:
            upregulated = add_upreg_info(upregulated, 'EGG_vs_Gp_MALE_Gp_EGG_UP')

        if gene in EGG_vs_Gp_MALE_Gp_MALE_UP:
            upregulated = add_upreg_info(upregulated, 'EGG_vs_Gp_MALE_Gp_MALE_UP')

        if gene in J2_vs_Gp_MALE_Gp_J2_UP:
            upregulated = add_upreg_info(upregulated, 'J2_vs_Gp_MALE_Gp_J2_UP')

        if gene in J2_vs_Gp_MALE_Gp_MALE_UP:
            upregulated = add_upreg_info(upregulated, 'J2_vs_Gp_MALE_Gp_MALE_UP')
        dup = gene_duplication[gene.rstrip()]
        #outfmt = "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene, motif_counts, 
                                                  # spry, Secreted, 
                                                   #annotation, GH, dup, up_J2_v_14dpi, up_J2_v_egg, upregulated)
        outfmt = "\t".join([gene, str(dogs), str(subs),
                            spry, Secreted, 
                            annotation, GH, dup, up_J2_v_14dpi, up_J2_v_egg,
                            Up_dpi14_vs_Gp_J2,
                            Up_dpi14_vs_Gp_EGG,
                            up_male_vs_dpi14,
                            upregulated])
        f_out.write(outfmt + "\n")
    f_out.close()


phobius = """GPALN_016360-T2	0	Y	n10-20c24/25o
GPALN_013463-T1	0	Y	n6-17c22/23o
GPALN_007409-T1	0	Y	n2-12c20/21o
GPALN_003911-T1	0	Y	n4-15c20/21o
GPALN_001146-T1	0	Y	n3-18c23/24o
GPALN_011857-T1	0	Y	n4-19c24/25o
GPALN_003912-T1	0	Y	n10-21c26/27o
GPALN_006913-T1	0	Y	n3-14c21/22o
GPALN_008092-T1	0	Y	n8-19c27/28o
GPALN_004683-T1	2	Y	n7-14c19/20o210-232i253-277o
GPALN_003832-T1	0	Y	n2-13c18/19o
GPALN_004168-T1	1	Y	n5-16c21/22o330-356i
GPALN_002690-T1	1	Y	n8-19c24/25o86-108i
GPALN_007200-T1	1	Y	n8-16c21/22o133-152i
GPALN_007170-T1	0	Y	n4-15c23/24o
GPALN_015947-T1	0	Y	n7-18c26/27o
GPALN_007350-T1	0	Y	n5-12c16/17o
GPALN_008515-T1	1	Y	n4-15c20/21o238-261i
GPALN_003622-T1	0	Y	n18-29c37/38o
GPALN_003313-T1	0	Y	n18-28c32/33o
GPALN_014976-T1	0	Y	n3-18c23/24o
GPALN_012627-T1	0	Y	n17-28c32/33o
GPALN_000702-T1	0	Y	n10-23c30/31o
GPALN_005665-T1	3	Y	n3-11c16/17o26-47i59-78o90-109i
GPALN_007597-T1	0	Y	n6-17c22/23o
GPALN_007923-T1	0	Y	n3-14c22/23o
GPALN_004487-T1	0	Y	n2-13c21/22o
GPALN_002218-T1	0	Y	n3-18c23/24o
GPALN_013563-T1	0	Y	n4-17c22/23o
GPALN_002356-T1	15	Y	n4-15c19/20o408-428i440-458o464-483i509-528o534-550i562-581o587-605i617-636o642-661i666-684o690-706i748-768o788-813i825-844o856-877i
GPALN_007201-T1	0	Y	n7-15c20/21o
GPALN_001093-T1	0	Y	n5-16c24/25o
GPALN_000982-T1	0	Y	n8-18c22/23o
GPALN_011605-T1	0	Y	n2-12c17/18o
GPALN_003969-T1	1	Y	n4-23c31/32o131-154i
GPALN_010990-T1	0	Y	n2-10c15/16o
GPALN_011561-T1	0	Y	n2-13c17/18o
GPALN_007139-T1	0	Y	n18-29c33/34o
GPALN_004023-T1	0	Y	n8-15c33/34o
GPALN_000928-T1	0	Y	n5-16c24/25o
GPALN_009556-T1	0	Y	n6-17c25/26o
GPALN_007058-T1	0	Y	n4-15c23/24o
GPALN_004019-T2	0	Y	n3-14c19/20o
GPALN_004665-T1	0	Y	n8-19c24/25o
GPALN_001405-T1	0	Y	n14-25c30/31o
GPALN_012150-T1	0	Y	n3-14c26/27o
GPALN_003855-T1	0	Y	n6-17c25/26o
GPALN_008910-T1	1	Y	n6-17c25/26o179-199i
GPALN_013623-T1	0	Y	n12-22c27/28o
GPALN_015738-T1	0	Y	n26-34c41/42o
GPALN_005309-T1	0	Y	n14-25c33/34o
GPALN_002782-T2	1	Y	n27-35c39/40o436-460i
GPALN_003578-T1	1	Y	n23-34c39/40o104-130i
GPALN_014941-T1	1	Y	n6-17c22/23o60-82i
GPALN_003641-T1	1	Y	n6-17c24/25o659-680i
GPALN_007871-T1	0	Y	n5-15c20/21o
GPALN_002152-T1	1	Y	n8-19c29/30o274-293i
GPALN_001641-T1	0	Y	n17-28c33/34o
GPALN_012821-T1	0	Y	n10-21c27/28o
GPALN_007710-T1	0	Y	n3-10c15/16o
GPALN_014464-T1	0	Y	n21-29c34/35o
GPALN_009695-T1	0	Y	n6-19c24/25o
GPALN_011651-T1	0	Y	n12-23c27/28o
GPALN_000323-T1	0	Y	n8-26c31/32o
GPALN_001215-T1	0	Y	n12-23c32/33o
GPALN_003771-T1	1	Y	n8-16c22/23o239-260i
GPALN_002702-T1	3	Y	n11-21c26/27o199-222i234-255o275-295i
GPALN_011996-T1	0	Y	n22-30c35/36o
GPALN_014576-T1	0	Y	n2-13c20/21o
GPALN_002510-T1	1	Y	n8-23c28/29o467-496i
GPALN_013263-T1	0	Y	n3-12c20/21o
GPALN_005239-T1	1	Y	n5-16c21/22o80-105i
GPALN_010127-T1	0	Y	n4-15c20/21o
GPALN_009094-T1	0	Y	n6-17c24/25o
GPALN_006755-T1	0	Y	n11-22c27/28o
GPALN_005209-T1	2	Y	n14-25c30/31o159-180i192-212o
GPALN_002370-T1	0	Y	n5-15c23/24o
GPALN_014857-T1	0	Y	n4-16c20/21o
GPALN_007681-T1	1	Y	n6-13c23/24o61-80i
GPALN_012973-T1	1	Y	n2-9c14/15o2113-2133i
GPALN_007238-T1	0	Y	n4-12c17/18o
GPALN_000440-T2	1	Y	n6-14c19/20o353-371i
GPALN_003741-T1	0	Y	n7-18c24/25o
GPALN_009909-T1	0	Y	n7-17c22/23o
GPALN_001887-T1	7	Y	n5-16c22/23o85-113i125-150o156-184i205-226o246-267i408-428o440-461i
GPALN_006042-T1	0	Y	n5-16c21/22o
GPALN_010159-T1	0	Y	n3-12c20/21o
GPALN_016087-T2	1	Y	n12-23c28/29o236-255i
GPALN_002180-T1	0	Y	n6-17c21/22o
GPALN_002788-T1	1	Y	n10-21c29/30o506-529i
GPALN_007216-T1	0	Y	n6-17c22/23o
GPALN_001532-T1	0	Y	n3-14c19/20o
GPALN_004347-T1	0	Y	n2-13c20/21o
GPALN_005591-T1	0	Y	n5-16c21/22o
GPALN_012295-T1	1	Y	n4-12c20/21o137-157i
GPALN_003635-T1	8	Y	n8-18c26/27o57-76i88-110o122-141i206-226o246-268i321-339o345-365i421-446o
GPALN_008141-T1	0	Y	n3-14c19/20o
GPALN_011592-T1	5	Y	n2-10c18/19o37-55i76-100o134-158i179-198o210-231i
GPALN_006124-T1	0	Y	n10-17c29/30o
GPALN_010727-T1	0	Y	n4-15c23/24o
GPALN_002151-T1	0	Y	n16-34c41/42o
GPALN_000071-T1	13	Y	n10-20c28/29o85-104i116-134o146-163i175-197o209-230i239-261o281-306i318-334o340-360i372-394o400-419i439-456o468-492i
GPALN_002602-T1	0	Y	n4-15c20/21o
GPALN_000307-T1	0	Y	n10-21c29/30o
GPALN_007673-T1	0	Y	n6-17c22/23o
GPALN_015245-T1	0	Y	n6-21c25/26o
GPALN_002595-T1	0	Y	n3-12c17/18o
GPALN_004399-T1	0	Y	n2-12c20/21o
GPALN_008150-T1	1	Y	n5-15c20/21o720-736i
GPALN_014870-T1	2	Y	n10-20c25/26o257-279i291-311o
GPALN_007543-T1	0	Y	n4-14c19/20o
GPALN_008978-T1	0	Y	n7-18c22/23o
GPALN_008471-T1	0	Y	n2-9c14/15o
GPALN_002181-T1	0	Y	n3-22c31/32o
GPALN_015291-T1	0	Y	n5-16c21/22o
GPALN_015697-T1	0	Y	n8-18c23/24o
GPALN_006046-T1	2	Y	n9-21c26/27o80-100i112-132o
GPALN_013667-T1	1	Y	n23-33c42/43o52-74i
GPALN_006057-T1	0	Y	n11-21c29/30o
GPALN_003838-T1	0	Y	n8-19c24/25o
GPALN_012986-T1	0	Y	n3-13c18/19o
GPALN_014697-T2	0	Y	n3-13c18/19o
GPALN_014471-T1	5	Y	n11-22c27/28o90-108i115-136o176-197i273-291o311-329i
GPALN_010120-T1	0	Y	n4-14c19/20o
GPALN_004972-T1	0	Y	n6-17c21/22o
GPALN_003083-T1	1	Y	n13-23c28/29o143-163i
GPALN_009899-T1	0	Y	n5-16c23/24o
GPALN_011879-T1	0	Y	n8-19c27/28o
GPALN_009017-T1	0	Y	n43-54c63/64o
GPALN_005443-T1	0	Y	n4-14c19/20o
GPALN_000385-T1	1	Y	n7-15c19/20o336-359i
GPALN_008251-T1	0	Y	n7-17c35/36o
GPALN_014529-T1	0	Y	n10-20c28/29o
GPALN_003907-T2	0	Y	n4-17c22/23o
GPALN_005016-T1	0	Y	n5-17c22/23o
GPALN_003903-T1	0	Y	n5-13c18/19o
GPALN_013136-T1	0	Y	n4-16c21/22o
GPALN_016191-T1	0	Y	n3-14c22/23o
GPALN_013448-T1	0	Y	n4-15c20/21o
GPALN_014665-T1	0	Y	n18-29c37/38o
GPALN_016257-T1	0	Y	n8-19c37/38o
GPALN_008086-T2	0	Y	n3-11c19/20o
GPALN_011506-T1	0	Y	n9-17c21/22o
GPALN_002562-T1	0	Y	n10-20c25/26o
GPALN_005903-T1	0	Y	n6-16c25/26o
GPALN_016172-T1	1	Y	n6-21c29/30o82-103i
GPALN_001041-T1	1	Y	n4-12c16/17o210-231i
GPALN_010517-T1	0	Y	n4-15c20/21o
GPALN_014238-T1	0	Y	n2-10c18/19o
GPALN_003955-T1	0	Y	n4-15c20/21o
GPALN_002448-T1	0	Y	n3-18c25/26o
GPALN_004926-T1	0	Y	n6-18c22/23o
GPALN_011240-T1	0	Y	n8-19c27/28o
GPALN_011894-T1	0	Y	n4-15c23/24o
GPALN_004880-T1	0	Y	n8-18c25/26o
GPALN_004045-T1	1	Y	n9-20c24/25o154-171i
GPALN_007795-T1	0	Y	n3-18c23/24o
GPALN_002252-T1	1	Y	n4-15c21/22o389-410i
GPALN_007180-T1	5	Y	n2-13c21/22o116-139i151-176o188-210i241-261o267-287i
GPALN_014698-T1	0	Y	n4-14c21/22o
GPALN_004700-T1	0	Y	n9-19c24/25o
GPALN_003416-T1	0	Y	n5-13c18/19o
GPALN_001208-T1	0	Y	n3-11c29/30o
GPALN_001482-T1	1	Y	n13-26c31/32o319-343i
GPALN_016380-T1	0	Y	n3-18c23/24o
GPALN_006067-T1	0	Y	n14-24c32/33o
GPALN_004020-T1	0	Y	n6-17c21/22o
GPALN_011260-T1	0	Y	n7-17c25/26o
GPALN_011957-T1	0	Y	n11-18c25/26o
GPALN_002043-T1	0	Y	n14-29c37/38o
GPALN_015833-T1	1	Y	n7-16c20/21o576-594i
GPALN_000168-T1	0	Y	n42-52c64/65o
GPALN_002044-T1	0	Y	n5-12c17/18o
GPALN_000157-T1	0	Y	n3-12c17/18o
GPALN_015281-T1	0	Y	n3-14c19/20o
GPALN_000311-T1	0	Y	n6-17c21/22o
GPALN_000490-T1	0	Y	n8-19c24/25o
GPALN_013242-T1	0	Y	n7-18c23/24o
GPALN_009896-T1	0	Y	n2-12c17/18o
GPALN_015248-T1	0	Y	n2-17c22/23o
GPALN_016042-T1	0	Y	n3-14c18/19o
GPALN_003971-T1	0	Y	n6-14c22/23o
GPALN_004396-T1	5	Y	n2-17c22/23o32-50i71-89o95-113i133-155o161-183i
GPALN_001016-T2	0	Y	n20-30c35/36o
GPALN_005068-T1	0	Y	n6-17c21/22o
GPALN_016107-T1	1	Y	n17-27c38/39o814-837i
GPALN_014371-T1	0	Y	n5-16c21/22o
GPALN_006591-T1	1	Y	n5-16c24/25o40-66i
GPALN_007620-T1	0	Y	n6-16c21/22o
GPALN_005752-T1	0	Y	n5-16c21/22o
GPALN_012706-T1	0	Y	n11-22c30/31o
GPALN_005577-T1	1	Y	n3-18c26/27o113-132i
GPALN_001062-T1	0	Y	n5-11c16/17o
GPALN_011617-T1	0	Y	n6-19c25/26o
GPALN_001042-T1	0	Y	n3-11c15/16o
GPALN_007177-T1	1	Y	n2-10c18/19o221-243i
GPALN_011535-T1	1	Y	n3-18c23/24o143-161i
GPALN_013366-T1	0	Y	n6-17c21/22o
GPALN_003683-T1	0	Y	n4-15c20/21o
GPALN_014237-T1	0	Y	n7-15c20/21o
GPALN_000397-T1	0	Y	n5-15c26/27o
GPALN_003233-T1	0	Y	n547-557c562/563o
GPALN_014706-T1	0	Y	n4-12c20/21o
GPALN_010671-T1	0	Y	n11-21c29/30o
GPALN_010628-T1	0	Y	n14-24c32/33o
GPALN_001912-T1	0	Y	n6-17c22/23o
GPALN_001026-T1	0	Y	n3-13c18/19o
GPALN_001625-T1	4	Y	n9-20c25/26o80-103i231-251o2098-2119i2140-2159o
GPALN_010457-T1	0	Y	n4-14c24/25o
GPALN_014931-T1	0	Y	n5-15c20/21o
GPALN_001286-T1	0	Y	n4-16c21/22o
GPALN_008431-T1	0	Y	n4-14c19/20o
GPALN_004667-T2	0	Y	n20-30c35/36o
GPALN_014919-T1	0	Y	n2-12c16/17o
GPALN_004977-T1	0	Y	n8-19c24/25o
GPALN_015831-T1	0	Y	n7-15c20/21o
GPALN_013264-T1	0	Y	n7-18c29/30o
GPALN_004463-T1	0	Y	n4-14c19/20o
GPALN_000485-T1	4	Y	n6-19c27/28o223-246i258-278o290-311i524-548o
GPALN_007443-T1	0	Y	n5-18c23/24o
GPALN_004593-T1	0	Y	n5-17c22/23o
GPALN_006914-T1	0	Y	n3-11c16/17o
GPALN_004251-T1	0	Y	n40-49c54/55o
GPALN_003055-T1	0	Y	n7-18c27/28o
GPALN_008951-T1	0	Y	n6-14c19/20o
GPALN_005713-T1	0	Y	n7-17c21/22o
GPALN_004506-T1	0	Y	n6-19c24/25o
GPALN_007639-T1	0	Y	n6-14c22/23o
GPALN_013382-T1	0	Y	n142-150c168/169o
GPALN_012297-T1	0	Y	n11-21c26/27o
GPALN_000111-T1	0	Y	n13-24c29/30o
GPALN_013214-T1	0	Y	n2-12c17/18o
GPALN_002651-T1	0	Y	n5-16c24/25o
GPALN_005007-T1	1	Y	n22-33c43/44o157-178i
GPALN_010805-T1	2	Y	n8-19c27/28o487-503i524-545o
GPALN_014408-T1	1	Y	n4-16c21/22o665-689i
GPALN_002508-T1	0	Y	n6-17c23/24o
GPALN_005024-T1	1	Y	n13-24c32/33o696-714i
GPALN_004481-T1	1	Y	n6-18c22/23o68-86i
GPALN_008317-T1	1	Y	n4-12c16/17o1443-1463i
GPALN_005912-T1	0	Y	n4-14c19/20o
GPALN_010289-T1	0	Y	n8-16c21/22o
GPALN_003867-T1	0	Y	n2-13c20/21o
GPALN_004955-T1	2	Y	n10-18c27/28o43-66i207-229o
GPALN_011618-T1	0	Y	n4-14c18/19o
GPALN_015423-T1	0	Y	n4-18c23/24o
GPALN_010659-T1	0	Y	n2-10c21/22o
GPALN_009154-T1	0	Y	n3-10c15/16o
GPALN_000536-T1	0	Y	n7-22c30/31o
GPALN_002433-T1	1	Y	n8-19c24/25o283-306i
GPALN_005528-T1	7	Y	n3-14c18/19o28-46i53-72o78-95i102-121o127-149i170-189o195-213i
GPALN_001982-T1	0	Y	n12-23c28/29o
GPALN_011891-T2	0	Y	n8-23c30/31o
GPALN_001849-T1	0	Y	n5-16c20/21o
GPALN_015769-T1	0	Y	n4-16c21/22o
GPALN_006433-T1	0	Y	n8-16c21/22o
GPALN_001938-T1	0	Y	n5-13c18/19o
GPALN_008079-T1	0	Y	n9-22c27/28o
GPALN_008457-T1	0	Y	n5-16c21/22o
GPALN_010910-T1	0	Y	n6-17c21/22o
GPALN_009323-T1	0	Y	n4-16c24/25o
GPALN_005014-T1	0	Y	n2-12c17/18o
GPALN_013984-T1	0	Y	n17-32c40/41o
GPALN_011411-T1	0	Y	n4-15c25/26o
GPALN_007702-T1	0	Y	n2-9c13/14o
GPALN_006283-T1	9	Y	n9-19c24/25o231-256i384-405o434-452i459-481o501-526i538-556o576-601i613-634o646-663i
GPALN_003488-T1	0	Y	n9-17c22/23o
GPALN_007564-T1	0	Y	n6-17c22/23o
GPALN_004456-T1	0	Y	n5-15c22/23o
GPALN_010642-T1	0	Y	n14-25c29/30o
GPALN_008391-T1	0	Y	n22-33c38/39o
GPALN_003772-T1	0	Y	n3-14c18/19o
GPALN_004119-T1	0	Y	n13-21c26/27o
GPALN_015952-T1	0	Y	n10-18c36/37o
GPALN_000668-T1	0	Y	n10-21c26/27o
GPALN_001258-T1	0	Y	n6-19c24/25o
GPALN_015649-T1	1	Y	n11-26c34/35o44-64i
GPALN_000253-T1	1	Y	n28-38c43/44o1142-1164i
GPALN_009490-T1	0	Y	n12-22c32/33o
GPALN_009721-T1	0	Y	n9-20c32/33o
GPALN_008154-T1	0	Y	n2-13c18/19o
GPALN_005531-T1	0	Y	n7-18c26/27o
GPALN_005433-T1	1	Y	n13-24c32/33o48-69i
GPALN_000136-T2	0	Y	n9-19c26/27o
GPALN_005864-T1	7	Y	n4-14c32/33o42-65i77-100o131-155i176-194o226-246i273-298o318-338i
GPALN_000434-T1	0	Y	n5-15c20/21o
GPALN_013142-T1	0	Y	n6-17c35/36o
GPALN_001505-T1	0	Y	n7-18c24/25o
GPALN_001182-T1	0	Y	n5-12c17/18o
GPALN_005866-T1	0	Y	n13-21c27/28o
GPALN_015184-T1	0	Y	n4-22c26/27o
GPALN_000569-T1	0	Y	n8-19c27/28o
GPALN_007934-T1	8	Y	n2-10c17/18o455-473i671-693o752-770i930-952o958-980i992-1013o1025-1041i1251-1274o
GPALN_009020-T1	6	Y	n16-27c32/33o48-75i111-130o142-164i258-277o283-307i328-354o
GPALN_001203-T1	0	Y	n4-14c19/20o
GPALN_002248-T1	0	Y	n12-23c28/29o
GPALN_007633-T1	0	Y	n3-14c19/20o
GPALN_005138-T1	0	Y	n3-11c16/17o
GPALN_003437-T1	13	Y	n4-15c33/34o139-162i262-288o300-319i331-353o359-382i394-416o436-457i492-512o532-556i568-588o624-652i664-685o691-711i
GPALN_012455-T2	0	Y	n3-21c25/26o
GPALN_008042-T1	0	Y	n7-19c23/24o
GPALN_014745-T1	0	Y	n6-14c22/23o
GPALN_009625-T1	2	Y	n6-17c29/30o754-773i794-812o
GPALN_011714-T1	0	Y	n5-16c21/22o
GPALN_004758-T1	2	Y	n3-14c18/19o28-49i61-87o
GPALN_009377-T2	1	Y	n10-21c39/40o549-571i
GPALN_006112-T1	0	Y	n6-19c24/25o
GPALN_009177-T1	0	Y	n13-24c31/32o
GPALN_010463-T1	0	Y	n4-15c26/27o
GPALN_008101-T1	0	Y	n3-13c31/32o
GPALN_003779-T1	0	Y	n2-13c31/32o
GPALN_008006-T1	2	Y	n6-19c23/24o1070-1089i1272-1297o
GPALN_010829-T1	0	Y	n4-19c28/29o
GPALN_002179-T1	0	Y	n4-16c25/26o
GPALN_002880-T2	0	Y	n4-15c23/24o
GPALN_010655-T1	1	Y	n2-13c19/20o29-45i
GPALN_003436-T1	0	Y	n14-22c33/34o
GPALN_009477-T1	0	Y	n5-15c19/20o
GPALN_010602-T1	0	Y	n2-13c18/19o
GPALN_003712-T1	0	Y	n8-19c24/25o
GPALN_015819-T1	3	Y	n5-13c21/22o91-114i126-149o169-198i
GPALN_015696-T1	1	Y	n10-21c28/29o205-226i
GPALN_001428-T1	9	Y	n4-14c19/20o29-50i71-91o97-118i125-144o150-169i189-208o228-251i258-278o284-303i
GPALN_005214-T1	9	Y	n5-16c21/22o349-370i382-402o544-566i573-600o1045-1065i1104-1127o1175-1192i1204-1228o1240-1256i
GPALN_008305-T1	2	Y	n3-14c19/20o38-67i87-106o
GPALN_014591-T1	0	Y	n2-12c20/21o
GPALN_014926-T1	0	Y	n23-33c41/42o
GPALN_002379-T1	0	Y	n5-16c21/22o
GPALN_002532-T1	0	Y	n3-14c22/23o
GPALN_008664-T1	0	Y	n8-19c24/25o
GPALN_003371-T1	0	Y	n12-24c28/29o
GPALN_005129-T1	0	Y	n17-28c36/37o
GPALN_006039-T2	0	Y	n11-18c23/24o
GPALN_010907-T1	0	Y	n3-14c19/20o
GPALN_000538-T1	0	Y	n11-22c30/31o
GPALN_003470-T1	0	Y	n3-14c19/20o
GPALN_016161-T1	0	Y	n6-19c28/29o
GPALN_014355-T1	0	Y	n10-21c25/26o
GPALN_005812-T1	0	Y	n5-14c19/20o
GPALN_009500-T1	0	Y	n3-18c23/24o
GPALN_002295-T1	0	Y	n5-20c24/25o
GPALN_009834-T1	0	Y	n39-50c55/56o
GPALN_014874-T1	0	Y	n6-17c22/23o
GPALN_004073-T1	0	Y	n5-18c22/23o
GPALN_004151-T1	0	Y	n4-13c17/18o
GPALN_015183-T1	0	Y	n8-18c26/27o
GPALN_012602-T1	0	Y	n15-25c30/31o
GPALN_006818-T1	0	Y	n5-20c24/25o
GPALN_013922-T1	0	Y	n4-11c16/17o
GPALN_005717-T1	0	Y	n6-16c24/25o
GPALN_006704-T1	0	Y	n4-14c18/19o
GPALN_008852-T1	2	Y	n6-16c28/29o558-578i585-605o
GPALN_009615-T1	1	Y	n8-19c31/32o644-669i
GPALN_013767-T1	3	Y	n4-15c20/21o164-183i195-216o751-770i
GPALN_004470-T1	0	Y	n5-16c21/22o
GPALN_009514-T1	1	Y	n6-19c24/25o48-69i
GPALN_002315-T1	1	Y	n8-19c24/25o123-149i
GPALN_013939-T1	0	Y	n6-18c23/24o
GPALN_010552-T1	0	Y	n4-17c22/23o
GPALN_012056-T1	0	Y	n4-15c20/21o
GPALN_004057-T1	0	Y	n9-20c24/25o
GPALN_015060-T1	0	Y	n13-21c27/28o
GPALN_015217-T1	4	Y	n10-18c24/25o250-268i275-296o308-329i534-553o
GPALN_007191-T1	0	Y	n4-13c18/19o
GPALN_004050-T1	0	Y	n4-12c17/18o
GPALN_007561-T1	1	Y	n7-17c24/25o268-289i
GPALN_005501-T1	0	Y	n4-14c24/25o
GPALN_004934-T1	1	Y	n4-11c16/17o459-477i
GPALN_004843-T1	0	Y	n6-17c24/25o
GPALN_012843-T1	7	Y	n8-19c24/25o67-91i103-122o149-168i189-209o229-254i370-391o406-427i
GPALN_013546-T1	0	Y	n25-36c41/42o
GPALN_007377-T1	0	Y	n7-19c24/25o
GPALN_007604-T1	0	Y	n6-17c24/25o
GPALN_012629-T1	1	Y	n5-15c19/20o398-419i
GPALN_011202-T1	10	Y	n3-14c19/20o29-49i61-80o92-113i125-142o178-199i281-303o315-336i343-365o377-396i416-437o
GPALN_013347-T1	0	Y	n4-15c20/21o
GPALN_005311-T1	0	Y	n2-9c14/15o
GPALN_011738-T1	0	Y	n3-14c19/20o
GPALN_014626-T1	0	Y	n3-14c21/22o
GPALN_008955-T1	1	Y	n16-27c35/36o45-66i
GPALN_014042-T1	0	Y	n4-19c23/24o
GPALN_002521-T1	0	Y	n6-16c21/22o
GPALN_015298-T1	0	Y	n2-20c24/25o
GPALN_006782-T1	0	Y	n5-15c23/24o
GPALN_003440-T1	0	Y	n10-21c29/30o
GPALN_004601-T1	0	Y	n5-12c17/18o
GPALN_016116-T1	1	Y	n5-18c23/24o220-244i
GPALN_006154-T1	0	Y	n9-24c29/30o
GPALN_000153-T3	0	Y	n20-33c38/39o
GPALN_014889-T1	1	Y	n6-18c23/24o39-60i
GPALN_007172-T1	0	Y	n7-15c20/21o
GPALN_013946-T1	0	Y	n23-36c54/55o
GPALN_005773-T1	0	Y	n6-17c22/23o
GPALN_005211-T1	0	Y	n38-49c60/61o
GPALN_014308-T1	0	Y	n5-16c21/22o
GPALN_000375-T1	3	Y	n5-20c25/26o111-131i322-341o353-369i
GPALN_008462-T1	0	Y	n3-14c19/20o
GPALN_008820-T1	0	Y	n12-31c35/36o
GPALN_008511-T1	0	Y	n5-12c20/21o
GPALN_011069-T1	0	Y	n25-34c42/43o
GPALN_000048-T1	0	Y	n10-21c29/30o
GPALN_013581-T1	0	Y	n6-17c22/23o
GPALN_010350-T1	1	Y	n4-15c20/21o97-115i
GPALN_012828-T1	0	Y	n3-18c25/26o
GPALN_005432-T1	1	Y	n7-15c19/20o280-305i
GPALN_012727-T1	1	Y	n7-17c22/23o418-438i
GPALN_013039-T1	0	Y	n2-12c17/18o
GPALN_014583-T1	1	Y	n4-17c25/26o200-224i
GPALN_009836-T1	0	Y	n4-15c20/21o
GPALN_009796-T1	0	Y	n9-20c24/25o
GPALN_003141-T1	0	Y	n9-19c26/27o
GPALN_015119-T1	0	Y	n3-15c20/21o
GPALN_002902-T1	1	Y	n18-26c44/45o87-112i
GPALN_005904-T1	0	Y	n5-16c21/22o
GPALN_006795-T1	0	Y	n4-14c22/23o
GPALN_008967-T1	1	Y	n17-28c33/34o744-764i
GPALN_013571-T1	0	Y	n8-16c21/22o
GPALN_013163-T1	1	Y	n4-15c20/21o392-412i
GPALN_008130-T1	0	Y	n9-19c23/24o
GPALN_002601-T1	0	Y	n9-20c25/26o
GPALN_000632-T1	0	Y	n17-25c30/31o
GPALN_007554-T1	0	Y	n8-20c25/26o
GPALN_004553-T1	0	Y	n6-16c21/22o
GPALN_012802-T1	5	Y	n3-11c16/17o26-48i68-91o97-114i134-152o164-184i
GPALN_013514-T1	1	Y	n19-30c34/35o492-514i
GPALN_006047-T1	2	Y	n9-21c26/27o76-96i108-125o
GPALN_015999-T1	0	Y	n2-13c19/20o
GPALN_004568-T1	0	Y	n3-16c21/22o
GPALN_010013-T1	0	Y	n10-23c30/31o
GPALN_003496-T1	1	Y	n6-17c22/23o50-74i
GPALN_004797-T1	0	Y	n4-15c23/24o
GPALN_001589-T1	1	Y	n13-24c31/32o520-543i
GPALN_007578-T1	0	Y	n4-14c22/23o
GPALN_005645-T1	1	Y	n14-25c30/31o440-462i
GPALN_004599-T1	0	Y	n4-14c19/20o
GPALN_013578-T1	0	Y	n4-15c23/24o
GPALN_006756-T1	0	Y	n11-22c27/28o
GPALN_013781-T1	0	Y	n6-17c23/24o
GPALN_008100-T1	0	Y	n3-18c23/24o
GPALN_000953-T1	5	Y	n5-17c22/23o103-125i364-383o395-413i420-439o521-540i
GPALN_014189-T1	0	Y	n7-16c20/21o
GPALN_000541-T1	0	Y	n11-19c27/28o
GPALN_001879-T1	0	Y	n3-13c19/20o
GPALN_013755-T1	5	Y	n6-16c34/35o44-73i85-102o122-141i148-167o173-192i
GPALN_014232-T1	1	Y	n5-16c21/22o203-226i
GPALN_013550-T1	0	Y	n4-22c27/28o
GPALN_012556-T1	0	Y	n10-20c25/26o
GPALN_002290-T1	0	Y	n3-13c17/18o
GPALN_015258-T1	0	Y	n15-26c31/32o
GPALN_011724-T1	0	Y	n3-14c19/20o
GPALN_014095-T1	0	Y	n11-21c26/27o
GPALN_016035-T1	0	Y	n2-12c24/25o
GPALN_016182-T1	0	Y	n4-12c17/18o
GPALN_005419-T1	0	Y	n10-21c26/27o
GPALN_014299-T1	0	Y	n3-13c18/19o
GPALN_010121-T1	0	Y	n7-18c22/23o
GPALN_002481-T1	0	Y	n4-14c18/19o
GPALN_000123-T1	0	Y	n3-11c17/18o
GPALN_003370-T1	0	Y	n7-18c25/26o
GPALN_003612-T1	1	Y	n14-26c31/32o902-924i
GPALN_010685-T1	0	Y	n2-13c20/21o
GPALN_013387-T1	0	Y	n3-13c21/22o
GPALN_013545-T1	0	Y	n9-20c38/39o
GPALN_012603-T1	5	Y	n32-42c50/51o84-105i146-167o179-200i212-235o241-261i
GPALN_003566-T1	0	Y	n8-19c27/28o
GPALN_010895-T1	0	Y	n4-15c20/21o
GPALN_009877-T1	0	Y	n3-11c16/17o
GPALN_007163-T1	0	Y	n4-15c23/24o
GPALN_010140-T1	0	Y	n6-16c22/23o
GPALN_010620-T1	0	Y	n3-11c16/17o
GPALN_010899-T2	0	Y	n10-20c24/25o
GPALN_006856-T1	0	Y	n5-20c24/25o
GPALN_010621-T1	1	Y	n3-13c18/19o126-146i
GPALN_013104-T1	0	Y	n3-11c16/17o
GPALN_007552-T1	0	Y	n5-16c22/23o
GPALN_010526-T1	0	Y	n3-14c18/19o
GPALN_001970-T1	0	Y	n7-16c28/29o
GPALN_010348-T1	0	Y	n4-14c18/19o
GPALN_011989-T1	0	Y	n9-20c27/28o
GPALN_002182-T1	0	Y	n8-19c28/29o
GPALN_006522-T1	0	Y	n20-30c48/49o
GPALN_005091-T1	1	Y	n9-18c22/23o496-515i
GPALN_011615-T2	0	Y	n4-14c18/19o
GPALN_005970-T1	0	Y	n8-18c23/24o
GPALN_012641-T1	0	Y	n4-15c19/20o
GPALN_008107-T1	1	Y	n4-17c22/23o2863-2885i
GPALN_000455-T1	0	Y	n7-18c22/23o
GPALN_010664-T1	1	Y	n19-29c35/36o369-394i
GPALN_005744-T1	2	Y	n27-38c42/43o360-376i397-417o
GPALN_000040-T1	8	Y	n13-24c29/30o245-268i280-298o310-335i595-615o921-944i956-974o986-1009i1278-1301o
GPALN_006864-T1	0	Y	n9-20c24/25o
GPALN_009721-T2	0	Y	n9-20c32/33o
GPALN_015294-T1	0	Y	n5-16c24/25o
GPALN_004239-T2	0	Y	n8-16c21/22o
GPALN_016014-T1	0	Y	n7-16c21/22o
GPALN_006806-T1	0	Y	n2-13c17/18o
GPALN_007179-T1	1	Y	n2-10c18/19o238-262i
GPALN_015419-T1	0	Y	n7-16c20/21o
GPALN_005802-T1	0	Y	n4-15c19/20o
GPALN_004624-T1	0	Y	n7-18c22/23o
GPALN_007719-T1	0	Y	n8-18c23/24o
GPALN_011543-T1	0	Y	n3-10c15/16o
GPALN_000343-T1	0	Y	n3-14c19/20o
GPALN_000199-T1	0	Y	n8-18c29/30o
GPALN_005609-T1	0	Y	n9-18c23/24o
GPALN_004223-T1	0	Y	n12-23c31/32o
GPALN_013247-T1	0	Y	n7-18c23/24o
GPALN_006860-T1	0	Y	n4-15c24/25o
GPALN_013531-T1	0	Y	n18-29c37/38o
GPALN_009258-T1	2	Y	n5-13c23/24o434-452i657-680o
GPALN_012005-T1	0	Y	n4-15c33/34o
GPALN_015179-T1	0	Y	n10-21c26/27o
GPALN_010733-T1	0	Y	n3-13c17/18o
GPALN_011778-T1	0	Y	n9-23c28/29o
GPALN_000543-T1	1	Y	n4-19c26/27o1025-1049i
GPALN_007648-T1	0	Y	n7-18c27/28o
GPALN_008880-T1	0	Y	n4-15c23/24o
GPALN_000079-T1	1	Y	n5-16c21/22o2918-2939i
GPALN_006627-T1	0	Y	n8-19c24/25o
GPALN_000499-T1	4	Y	n18-26c35/36o1362-1379i1386-1405o1417-1439i1448-1471o
GPALN_008027-T1	1	Y	n2-10c22/23o377-399i
GPALN_004059-T1	1	Y	n9-20c26/27o433-455i
GPALN_012414-T1	0	Y	n3-14c19/20o
GPALN_000858-T2	2	Y	n8-16c28/29o395-421i610-629o
GPALN_007916-T1	7	Y	n5-16c34/35o172-197i209-234o246-268i289-308o339-365i686-704o724-744i
GPALN_008601-T1	0	Y	n11-22c27/28o
GPALN_013557-T1	0	Y	n3-10c15/16o
GPALN_011613-T1	0	Y	n6-19c25/26o
GPALN_001828-T1	0	Y	n3-11c19/20o
GPALN_003326-T1	0	Y	n3-18c23/24o
GPALN_006413-T1	0	Y	n15-26c37/38o
GPALN_014006-T1	2	Y	n3-15c19/20o236-253i274-293o
GPALN_013224-T1	0	Y	n2-12c17/18o
GPALN_010437-T1	0	Y	n7-16c21/22o
GPALN_016213-T1	1	Y	n7-18c23/24o496-516i
GPALN_009217-T1	0	Y	n9-19c24/25o
GPALN_010538-T1	0	Y	n8-19c31/32o
GPALN_006660-T1	1	Y	n4-15c20/21o330-347i
GPALN_015499-T2	1	Y	n9-20c29/30o373-393i
GPALN_011626-T1	0	Y	n5-16c25/26o
GPALN_004675-T1	0	Y	n5-17c22/23o
GPALN_005637-T1	11	Y	n3-12c16/17o405-425i498-521o527-547i559-583o589-608i674-691o697-715i722-743o749-777i789-812o818-839i
GPALN_015532-T1	0	Y	n7-21c26/27o
GPALN_007917-T1	0	Y	n8-15c20/21o
GPALN_010861-T1	0	Y	n3-11c16/17o
GPALN_003481-T1	0	Y	n12-20c25/26o
GPALN_013256-T1	1	Y	n6-17c22/23o160-186i
GPALN_005344-T1	0	Y	n4-15c20/21o
GPALN_001368-T1	9	Y	n4-14c18/19o222-242i285-310o322-344i356-380o386-409i443-465o477-503i515-535o547-569i
GPALN_002774-T1	0	Y	n20-31c36/37o
GPALN_014737-T1	0	Y	n11-19c24/25o
GPALN_009412-T1	0	Y	n12-23c29/30o
GPALN_013150-T1	0	Y	n9-16c24/25o
GPALN_014567-T1	0	Y	n7-15c20/21o
GPALN_015231-T1	3	Y	n3-11c16/17o62-84i125-142o162-181i
GPALN_004116-T1	0	Y	n4-14c19/20o
GPALN_007586-T1	13	Y	n6-16c20/21o67-85i97-115o127-143i155-174o180-199i220-242o262-285i292-308o314-333i345-365o377-405i417-435o447-471i
GPALN_001972-T1	1	Y	n9-21c27/28o345-368i
GPALN_003442-T1	0	Y	n2-13c18/19o
GPALN_000109-T1	0	Y	n8-21c26/27o
GPALN_005927-T1	0	Y	n8-18c23/24o
GPALN_013709-T1	4	Y	n10-19c24/25o34-52i73-90o116-136i148-169o
GPALN_012955-T1	0	Y	n6-17c22/23o
GPALN_009648-T1	1	Y	n17-27c34/35o1189-1211i
GPALN_008584-T1	9	Y	n5-16c22/23o227-247i336-356o368-388i467-485o491-507i519-539o608-630i637-662o668-687i
GPALN_006590-T1	0	Y	n3-14c18/19o
GPALN_003057-T1	0	Y	n8-19c24/25o
GPALN_010549-T1	0	Y	n7-18c23/24o
GPALN_003382-T1	1	Y	n8-19c24/25o139-161i
GPALN_014750-T1	0	Y	n2-10c18/19o
GPALN_011437-T1	0	Y	n6-17c25/26o
GPALN_003081-T1	0	Y	n7-14c19/20o
GPALN_011378-T1	0	Y	n7-17c22/23o
GPALN_005159-T1	0	Y	n3-14c20/21o
GPALN_000608-T1	1	Y	n4-15c21/22o910-936i
GPALN_004671-T1	0	Y	n4-15c20/21o
GPALN_001729-T1	0	Y	n8-26c31/32o
GPALN_013535-T1	0	Y	n3-14c19/20o
GPALN_003697-T2	1	Y	n2-12c17/18o33-55i
GPALN_003398-T1	0	Y	n12-23c27/28o
GPALN_006998-T1	1	Y	n11-22c30/31o124-142i
GPALN_006085-T1	0	Y	n7-17c22/23o
GPALN_002307-T1	0	Y	n2-13c18/19o
GPALN_008213-T1	1	Y	n3-10c14/15o929-947i
GPALN_000312-T1	4	Y	n22-32c40/41o628-646i658-675o681-700i776-797o
GPALN_001674-T1	0	Y	n38-51c59/60o
GPALN_009835-T1	0	Y	n2-13c20/21o
GPALN_014151-T1	0	Y	n12-22c27/28o
GPALN_012357-T1	0	Y	n3-18c23/24o
GPALN_007049-T1	0	Y	n4-15c20/21o
GPALN_005196-T1	0	Y	n5-16c20/21o
GPALN_011456-T1	0	Y	n2-13c19/20o
GPALN_010511-T1	0	Y	n3-14c18/19o
GPALN_002028-T1	1	Y	n11-19c28/29o836-860i
GPALN_000162-T1	0	Y	n4-15c20/21o
GPALN_015675-T1	0	Y	n8-18c23/24o
GPALN_005064-T1	0	Y	n3-14c19/20o
GPALN_015354-T1	0	Y	n8-19c25/26o
GPALN_004130-T1	0	Y	n7-18c25/26o
GPALN_011222-T1	1	Y	n36-49c56/57o378-397i
GPALN_014964-T1	0	Y	n6-16c24/25o
GPALN_010639-T1	0	Y	n8-19c26/27o
GPALN_005017-T1	3	Y	n6-13c17/18o567-588i678-701o876-898i
GPALN_008002-T1	0	Y	n9-20c27/28o
GPALN_016121-T1	1	Y	n8-16c21/22o715-736i
GPALN_008935-T1	0	Y	n8-16c21/22o
GPALN_012513-T1	0	Y	n7-18c27/28o
GPALN_000181-T1	0	Y	n6-17c22/23o
GPALN_015309-T1	0	Y	n5-20c24/25o
GPALN_001758-T1	0	Y	n7-18c22/23o
GPALN_002503-T1	13	Y	n2-10c15/16o31-47i68-87o107-131i210-230o300-324i400-421o433-454i466-488o500-521i607-632o675-695i707-729o735-753i
GPALN_012774-T1	0	Y	n9-20c28/29o
GPALN_001251-T1	0	Y	n2-13c18/19o
GPALN_002346-T1	2	Y	n6-17c22/23o778-795i807-829o
GPALN_013326-T1	0	Y	n4-15c20/21o
GPALN_006652-T1	0	Y	n3-15c20/21o
GPALN_004288-T1	0	Y	n3-16c24/25o
GPALN_014438-T1	0	Y	n5-17c21/22o
GPALN_003454-T1	0	Y	n4-16c28/29o
GPALN_006623-T1	0	Y	n4-11c19/20o
GPALN_007720-T1	1	Y	n8-18c22/23o115-139i
GPALN_005063-T1	0	Y	n6-16c24/25o
GPALN_009904-T1	0	Y	n7-17c22/23o
GPALN_005493-T1	5	Y	n26-33c37/38o47-68i80-97o103-123i144-162o168-185i
GPALN_005008-T1	0	Y	n4-15c27/28o
GPALN_003759-T1	1	Y	n11-18c26/27o123-141i
GPALN_014759-T1	0	Y	n2-13c21/22o
GPALN_003643-T1	0	Y	n9-18c23/24o
GPALN_005902-T1	0	Y	n5-16c22/23o
GPALN_002965-T1	0	Y	n3-13c17/18o
GPALN_007525-T1	1	Y	n3-11c16/17o157-174i
GPALN_014370-T1	0	Y	n5-16c21/22o
GPALN_012523-T1	0	Y	n8-19c24/25o
GPALN_009610-T1	1	Y	n4-15c21/22o277-294i
GPALN_015505-T1	0	Y	n7-17c25/26o
GPALN_000767-T1	0	Y	n15-26c32/33o
GPALN_007548-T1	0	Y	n3-14c19/20o
GPALN_010880-T1	0	Y	n6-17c25/26o
GPALN_005754-T1	0	Y	n7-18c23/24o
GPALN_004460-T1	0	Y	n7-18c23/24o
GPALN_003445-T1	1	Y	n6-17c22/23o441-464i
GPALN_006828-T1	0	Y	n3-13c17/18o
GPALN_006043-T1	2	Y	n9-21c26/27o69-91i103-126o
GPALN_015300-T1	0	Y	n8-18c25/26o
GPALN_014397-T1	0	Y	n3-13c17/18o
GPALN_006088-T2	0	Y	n3-15c23/24o
GPALN_010433-T1	0	Y	n7-15c20/21o
GPALN_008903-T1	0	Y	n16-27c32/33o
GPALN_014403-T1	0	Y	n4-16c21/22o
GPALN_005270-T1	0	Y	n4-15c21/22o
GPALN_015690-T1	0	Y	n9-17c22/23o
GPALN_009137-T1	1	Y	n10-21c25/26o739-758i
GPALN_007286-T1	0	Y	n5-16c20/21o
GPALN_005276-T1	0	Y	n7-18c23/24o
GPALN_012284-T1	0	Y	n8-19c24/25o
GPALN_011013-T1	0	Y	n7-25c29/30o
GPALN_006646-T1	0	Y	n4-15c23/24o
GPALN_013257-T1	0	Y	n9-17c25/26o
GPALN_007685-T1	0	Y	n7-22c27/28o
GPALN_002129-T2	7	Y	n3-13c21/22o451-474i494-511o531-560i673-694o734-755i775-794o800-821i
GPALN_008356-T1	3	Y	n7-17c21/22o81-102i123-143o163-189i
GPALN_007603-T1	1	Y	n11-26c31/32o233-255i
GPALN_004818-T1	0	Y	n7-18c26/27o
GPALN_007267-T1	0	Y	n5-15c23/24o
GPALN_007039-T1	0	Y	n7-18c22/23o
GPALN_010710-T1	0	Y	n7-18c23/24o
GPALN_012294-T1	0	Y	n5-16c22/23o
GPALN_011335-T1	7	Y	n10-21c26/27o238-259i271-294o300-319i331-357o377-397i409-431o443-466i
GPALN_009507-T1	0	Y	n10-21c28/29o
GPALN_005777-T1	0	Y	n3-16c21/22o
GPALN_013662-T1	0	Y	n2-11c16/17o
GPALN_003913-T1	0	Y	n10-21c26/27o
GPALN_012064-T1	0	Y	n5-20c24/25o
GPALN_002183-T1	0	Y	n4-15c20/21o
GPALN_007859-T1	1	Y	n11-21c26/27o662-685i
GPALN_001055-T1	0	Y	n17-28c33/34o
GPALN_002784-T1	0	Y	n19-30c35/36o
GPALN_002941-T1	1	Y	n13-24c29/30o218-242i
GPALN_001224-T1	0	Y	n18-29c33/34o
GPALN_006338-T1	0	Y	n9-18c26/27o
GPALN_001494-T1	1	Y	n5-16c21/22o114-132i
GPALN_000515-T1	1	Y	n21-30c48/49o86-104i
GPALN_012634-T2	0	Y	n34-45c52/53o
GPALN_000063-T1	1	Y	n4-11c15/16o727-749i
GPALN_012618-T1	0	Y	n6-14c22/23o
GPALN_007199-T1	0	Y	n4-14c22/23o
GPALN_009402-T1	0	Y	n5-16c21/22o
GPALN_014395-T1	0	Y	n4-15c24/25o
GPALN_004924-T1	0	Y	n6-16c20/21o
GPALN_011731-T1	0	Y	n28-39c44/45o
GPALN_002865-T1	7	Y	n17-28c33/34o118-143i155-177o197-219i231-253o293-317i338-365o385-410i
GPALN_008149-T1	0	Y	n14-25c33/34o
GPALN_004105-T1	0	Y	n20-28c33/34o
GPALN_003175-T1	1	Y	n4-15c23/24o226-247i
GPALN_008078-T1	0	Y	n14-27c32/33o
GPALN_015355-T1	0	Y	n7-18c24/25o
GPALN_006725-T1	1	Y	n6-17c22/23o184-207i
GPALN_003950-T1	0	Y	n4-15c20/21o
GPALN_005071-T1	0	Y	n13-24c33/34o
GPALN_007574-T2	0	Y	n3-11c18/19o
GPALN_010059-T1	0	Y	n7-17c25/26o
GPALN_002917-T1	0	Y	n8-27c34/35o
GPALN_000080-T1	0	Y	n14-25c30/31o
GPALN_002729-T1	1	Y	n11-18c23/24o856-884i
GPALN_007539-T1	7	Y	n9-20c25/26o620-642i663-681o693-713i734-762o782-810i850-874o894-915i
GPALN_008137-T1	0	Y	n3-13c20/21o
GPALN_016269-T1	0	Y	n10-20c27/28o
GPALN_001924-T1	0	Y	n8-17c21/22o
GPALN_005953-T1	0	Y	n3-18c23/24o
GPALN_005742-T1	0	Y	n8-19c27/28o
GPALN_014279-T1	0	Y	n7-20c25/26o
GPALN_015347-T1	0	Y	n7-18c24/25o
GPALN_006261-T1	0	Y	n6-17c35/36o
GPALN_013707-T1	1	Y	n7-18c23/24o85-106i
GPALN_001503-T1	1	Y	n35-46c51/52o622-643i
GPALN_015940-T1	0	Y	n16-27c31/32o
GPALN_005058-T1	0	Y	n4-15c19/20o
GPALN_001433-T1	0	Y	n14-22c27/28o
GPALN_014019-T1	0	Y	n13-21c28/29o
GPALN_002198-T1	5	Y	n4-11c15/16o31-56i68-87o118-141i162-185o223-241i
GPALN_000568-T1	0	Y	n2-9c14/15o
GPALN_008540-T1	4	Y	n4-15c23/24o236-253i262-282o294-316i503-522o
GPALN_007840-T1	1	Y	n19-30c38/39o241-265i
GPALN_011927-T1	1	Y	n19-29c33/34o1558-1579i
GPALN_008132-T1	0	Y	n4-19c24/25o
GPALN_002578-T1	0	Y	n6-17c24/25o
GPALN_012293-T1	1	Y	n2-12c20/21o136-156i
GPALN_001008-T1	1	Y	n12-23c30/31o729-752i
GPALN_014333-T1	0	Y	n8-22c27/28o
GPALN_014416-T1	0	Y	n6-16c22/23o
GPALN_007114-T1	0	Y	n3-14c19/20o
GPALN_003876-T1	0	Y	n2-12c16/17o
GPALN_006031-T1	0	Y	n11-18c23/24o
GPALN_006534-T1	0	Y	n4-15c19/20o
GPALN_007715-T1	0	Y	n3-13c18/19o
GPALN_005166-T1	0	Y	n3-14c19/20o
GPALN_014738-T1	0	Y	n19-30c35/36o
GPALN_004921-T1	1	Y	n12-23c28/29o237-262i
GPALN_012924-T1	1	Y	n17-28c33/34o180-198i
GPALN_000341-T1	0	Y	n3-14c18/19o
GPALN_004107-T1	0	Y	n4-12c17/18o
GPALN_016201-T1	0	Y	n8-19c25/26o
GPALN_003433-T1	1	Y	n6-16c24/25o523-542i
GPALN_014666-T1	0	Y	n7-17c23/24o
GPALN_010096-T1	0	Y	n15-27c32/33o
GPALN_001018-T1	0	Y	n4-12c16/17o
GPALN_004552-T1	0	Y	n6-16c21/22o
GPALN_004458-T1	0	Y	n4-15c20/21o
GPALN_012520-T1	12	Y	n5-13c18/19o28-46i67-86o106-127i134-151o171-189i201-222o228-251i318-335o341-357i369-391o411-441i453-479o
GPALN_010962-T1	1	Y	n6-17c29/30o66-89i
GPALN_010930-T1	0	Y	n24-32c40/41o
GPALN_012099-T1	0	Y	n3-14c32/33o
GPALN_012331-T1	0	Y	n10-21c27/28o
GPALN_002958-T1	0	Y	n13-21c26/27o
GPALN_010088-T1	0	Y	n5-14c19/20o
GPALN_002862-T1	0	Y	n3-14c19/20o
GPALN_008342-T1	1	Y	n10-21c27/28o159-178i
GPALN_003943-T1	0	Y	n25-36c41/42o
GPALN_010627-T1	2	Y	n6-19c27/28o70-95i116-137o
GPALN_004384-T1	0	Y	n3-13c21/22o
GPALN_000167-T1	1	Y	n9-20c32/33o1488-1510i
GPALN_004757-T1	0	Y	n7-15c20/21o
GPALN_015279-T1	0	Y	n5-16c21/22o
GPALN_002695-T1	0	Y	n4-12c20/21o
GPALN_008999-T1	0	Y	n6-17c21/22o
GPALN_007371-T1	0	Y	n5-16c21/22o
GPALN_013876-T1	0	Y	n24-35c47/48o
GPALN_005554-T1	0	Y	n2-12c17/18o
GPALN_014844-T1	0	Y	n3-11c16/17o
GPALN_008145-T1	0	Y	n2-12c20/21o
GPALN_008098-T1	0	Y	n3-18c23/24o
GPALN_015891-T1	0	Y	n6-17c22/23o
GPALN_005736-T1	3	Y	n2-12c17/18o1112-1131i1143-1160o1172-1190i
GPALN_009293-T1	0	Y	n8-19c24/25o
GPALN_004564-T1	0	Y	n10-21c26/27o
GPALN_004616-T1	4	Y	n3-14c19/20o261-283i290-312o324-346i602-619o
GPALN_009156-T1	1	Y	n3-14c19/20o151-174i
GPALN_010572-T1	11	Y	n6-16c20/21o57-79i91-111o117-138i150-173o179-199i271-288o308-329i336-355o367-387i408-427o433-453i
GPALN_003329-T1	0	Y	n2-9c14/15o
GPALN_005199-T1	1	Y	n10-21c26/27o166-185i
GPALN_008116-T1	1	Y	n9-17c21/22o45-70i
GPALN_011036-T1	1	Y	n13-23c28/29o155-188i
GPALN_013536-T1	0	Y	n4-13c18/19o
GPALN_012062-T1	0	Y	n3-14c32/33o
GPALN_004681-T1	0	Y	n2-11c16/17o
GPALN_002129-T1	7	Y	n3-13c21/22o451-474i494-511o531-560i673-694o734-755i775-794o800-821i
GPALN_003047-T1	1	Y	n13-25c30/31o249-275i
GPALN_001874-T1	2	Y	n7-20c25/26o307-326i687-706o
GPALN_014150-T1	0	Y	n2-12c17/18o
GPALN_011124-T1	2	Y	n4-14c19/20o29-46i170-189o
GPALN_005741-T1	0	Y	n4-15c20/21o
GPALN_013128-T1	0	Y	n7-15c20/21o
GPALN_015006-T1	0	Y	n7-18c26/27o
GPALN_010800-T1	0	Y	n6-13c18/19o
GPALN_002203-T1	0	Y	n2-10c14/15o
GPALN_007507-T1	0	Y	n4-15c21/22o
GPALN_002455-T1	0	Y	n4-12c17/18o
GPALN_013552-T2	0	Y	n12-27c35/36o
GPALN_005103-T1	0	Y	n4-15c20/21o
GPALN_014785-T1	0	Y	n3-13c18/19o
GPALN_007120-T1	0	Y	n18-29c33/34o
GPALN_007198-T3	0	Y	n2-9c14/15o
GPALN_000394-T1	0	Y	n5-16c21/22o
GPALN_003462-T1	0	Y	n4-14c18/19o
GPALN_000424-T1	0	Y	n10-17c22/23o
GPALN_013118-T1	0	Y	n5-12c17/18o
GPALN_015745-T1	1	Y	n6-18c23/24o166-188i
GPALN_014430-T1	0	Y	n4-15c20/21o
GPALN_007647-T1	0	Y	n4-15c20/21o
GPALN_002300-T1	0	Y	n5-20c24/25o
GPALN_008747-T2	1	Y	n4-12c17/18o27-45i
GPALN_012762-T1	1	Y	n4-15c20/21o174-203i
GPALN_001729-T2	0	Y	n8-26c31/32o
GPALN_003891-T1	0	Y	n8-21c26/27o
GPALN_015094-T1	0	Y	n6-16c20/21o
GPALN_015323-T1	0	Y	n2-10c15/16o
GPALN_016202-T3	7	Y	n14-33c42/43o244-264i276-297o317-335i347-366o386-408i444-464o495-520i
GPALN_003283-T1	0	Y	n8-19c23/24o
GPALN_006304-T1	0	Y	n2-17c25/26o
GPALN_007220-T1	0	Y	n7-18c23/24o
GPALN_000453-T1	2	Y	n4-12c17/18o300-323i355-374o
GPALN_015116-T1	0	Y	n2-13c18/19o
GPALN_006358-T1	4	Y	n16-27c33/34o189-208i220-239o273-292i304-324o
GPALN_007178-T1	1	Y	n2-10c18/19o238-262i
GPALN_011156-T1	0	Y	n8-19c24/25o
GPALN_012580-T1	1	Y	n9-19c23/24o193-214i
GPALN_005673-T1	0	Y	n2-12c20/21o
GPALN_011531-T1	0	Y	n4-14c19/20o
GPALN_001599-T1	0	Y	n5-12c16/17o
GPALN_007218-T1	0	Y	n7-18c23/24o
GPALN_016359-T1	0	Y	n6-17c29/30o
GPALN_005599-T1	0	Y	n4-18c23/24o
GPALN_012234-T1	0	Y	n8-17c22/23o
GPALN_005910-T1	0	Y	n4-14c19/20o
GPALN_011160-T1	0	Y	n10-20c25/26o
GPALN_009571-T1	2	Y	n7-18c23/24o55-75i87-108o
GPALN_009029-T1	2	Y	n8-18c25/26o476-496i625-651o
GPALN_005202-T1	0	Y	n6-19c23/24o
GPALN_005012-T1	0	Y	n7-18c23/24o
GPALN_003773-T1	0	Y	n6-16c21/22o
GPALN_012113-T1	0	Y	n14-24c28/29o
GPALN_008097-T1	0	Y	n3-18c23/24o
GPALN_000922-T1	0	Y	n3-11c29/30o
GPALN_013831-T1	0	Y	n16-27c35/36o
GPALN_014369-T1	0	Y	n5-16c21/22o
GPALN_015676-T1	0	Y	n8-18c23/24o
GPALN_015749-T1	0	Y	n45-56c64/65o
GPALN_000640-T1	0	Y	n11-22c27/28o
GPALN_004026-T1	0	Y	n10-23c28/29o
GPALN_004799-T1	0	Y	n8-19c28/29o
GPALN_015332-T1	0	Y	n7-18c26/27o
GPALN_006440-T1	0	Y	n8-18c23/24o
GPALN_012803-T1	0	Y	n18-29c33/34o
GPALN_003434-T2	0	Y	n16-27c32/33o
GPALN_013771-T1	0	Y	n2-13c18/19o
GPALN_015312-T1	0	Y	n3-10c18/19o
GPALN_004901-T1	0	Y	n5-20c24/25o
GPALN_012406-T1	0	Y	n3-14c19/20o
GPALN_003185-T1	0	Y	n4-19c24/25o
GPALN_006086-T1	0	Y	n7-17c22/23o
GPALN_003977-T1	0	Y	n3-13c18/19o
GPALN_007820-T1	0	Y	n11-18c23/24o
GPALN_010119-T1	1	Y	n4-19c24/25o155-177i
GPALN_002001-T1	0	Y	n3-14c19/20o
GPALN_010634-T1	11	Y	n11-18c23/24o39-61i73-95o101-124i144-162o182-205i241-266o278-297i309-333o339-359i380-398o410-431i
GPALN_004807-T1	0	Y	n14-25c29/30o
GPALN_004734-T1	0	Y	n9-20c24/25o
GPALN_010742-T1	0	Y	n4-15c22/23o
GPALN_014415-T1	0	Y	n10-17c22/23o
GPALN_014514-T1	0	Y	n14-25c30/31o
GPALN_002349-T1	0	Y	n5-16c24/25o
GPALN_003990-T1	0	Y	n5-13c21/22o
GPALN_012550-T1	1	Y	n3-21c27/28o164-188i
GPALN_015186-T1	0	Y	n3-16c21/22o
GPALN_014799-T1	0	Y	n5-14c18/19o
GPALN_004369-T1	0	Y	n5-15c23/24o
GPALN_016087-T1	1	Y	n12-23c28/29o239-259i
GPALN_008082-T1	1	Y	n4-15c20/21o302-325i
GPALN_002227-T1	0	Y	n3-11c16/17o
GPALN_000898-T1	0	Y	n21-31c36/37o
GPALN_008976-T1	1	Y	n16-27c33/34o49-70i
GPALN_003844-T1	0	Y	n2-9c14/15o
GPALN_005246-T1	2	Y	n6-13c24/25o199-220i241-258o
GPALN_011170-T1	0	Y	n14-24c29/30o
GPALN_004008-T1	0	Y	n4-14c22/23o
GPALN_000124-T1	0	Y	n4-15c23/24o
GPALN_015579-T1	0	Y	n5-16c21/22o
GPALN_002318-T1	1	Y	n4-15c20/21o165-186i
GPALN_009639-T1	0	Y	n6-17c22/23o
GPALN_004851-T1	0	Y	n28-38c43/44o
GPALN_000377-T1	0	Y	n11-22c29/30o
GPALN_007706-T1	0	Y	n7-16c24/25o
GPALN_002618-T1	0	Y	n6-14c19/20o
GPALN_011440-T1	8	Y	n9-17c25/26o266-284i296-315o321-340i347-367o379-399i420-439o445-462i474-492o
GPALN_008152-T1	0	Y	n4-16c21/22o
GPALN_002126-T1	1	Y	n3-15c20/21o424-442i
GPALN_001748-T1	0	Y	n10-19c24/25o
GPALN_011607-T1	0	Y	n6-19c25/26o
GPALN_000384-T1	1	Y	n19-29c34/35o348-371i
GPALN_011006-T1	0	Y	n18-29c34/35o
GPALN_014500-T1	0	Y	n11-18c23/24o
GPALN_003552-T1	0	Y	n15-25c43/44o
GPALN_011754-T1	1	Y	n8-19c25/26o161-181i
GPALN_012127-T1	0	Y	n8-19c24/25o
GPALN_009497-T1	0	Y	n3-13c18/19o
GPALN_012476-T1	0	Y	n15-26c31/32o
GPALN_003829-T1	2	Y	n5-17c22/23o1044-1066i1078-1104o
GPALN_010554-T1	0	Y	n8-18c25/26o
GPALN_001759-T1	0	Y	n6-17c22/23o
GPALN_000562-T1	2	Y	n5-16c21/22o99-124i335-360o
GPALN_002702-T2	1	Y	n11-21c26/27o199-222i
GPALN_012846-T1	0	Y	n3-11c19/20o
GPALN_015743-T1	1	Y	n18-29c37/38o322-347i
GPALN_009636-T1	0	Y	n6-17c22/23o
GPALN_002516-T1	0	Y	n6-17c24/25o
GPALN_014070-T1	0	Y	n9-20c25/26o
GPALN_013106-T1	0	Y	n3-18c23/24o
GPALN_004239-T1	0	Y	n8-16c21/22o
GPALN_015930-T1	1	Y	n21-31c39/40o143-163i
GPALN_003192-T1	0	Y	n7-17c21/22o
GPALN_008824-T1	0	Y	n10-21c29/30o
GPALN_003517-T1	0	Y	n37-48c56/57o
GPALN_004244-T1	2	Y	n5-16c28/29o768-792i891-908o
GPALN_015280-T1	0	Y	n4-14c19/20o
GPALN_003831-T1	0	Y	n8-19c26/27o
GPALN_006486-T1	0	Y	n9-19c24/25o
GPALN_011436-T1	0	Y	n9-20c31/32o
GPALN_003835-T1	0	Y	n8-19c24/25o
GPALN_013092-T1	0	Y	n8-15c33/34o
GPALN_014851-T1	0	Y	n3-14c19/20o
GPALN_014292-T1	0	Y	n4-13c18/19o
GPALN_006578-T1	0	Y	n8-20c25/26o
GPALN_005821-T1	0	Y	n4-15c21/22o
GPALN_010771-T1	1	Y	n3-14c18/19o34-58i
GPALN_013350-T1	0	Y	n4-15c20/21o
GPALN_015813-T1	0	Y	n18-26c31/32o
GPALN_007410-T1	1	Y	n6-17c22/23o243-264i
GPALN_000430-T1	3	Y	n2-13c18/19o648-669i751-773o931-955i
GPALN_010908-T1	0	Y	n7-18c23/24o
GPALN_003380-T1	0	Y	n2-13c17/18o
GPALN_002071-T1	0	Y	n8-19c27/28o
GPALN_015178-T1	0	Y	n10-21c26/27o
GPALN_005815-T1	0	Y	n4-15c20/21o
GPALN_006944-T1	1	Y	n10-25c30/31o546-567i
GPALN_009364-T1	1	Y	n8-19c24/25o371-398i
GPALN_014320-T1	0	Y	n8-22c27/28o
GPALN_014730-T1	1	Y	n7-20c25/26o146-165i
GPALN_001955-T1	8	Y	n6-17c22/23o812-840i847-865o885-904i935-956o976-1002i1023-1046o1058-1079i1506-1525o
GPALN_004651-T1	7	Y	n3-14c19/20o113-134i146-167o211-237i258-278o284-304i325-347o367-390i
GPALN_006294-T1	0	Y	n5-13c25/26o
GPALN_009396-T1	1	Y	n16-27c34/35o436-456i
GPALN_014755-T1	0	Y	n2-10c15/16o
GPALN_014645-T1	3	Y	n17-28c33/34o43-63i75-95o115-137i
GPALN_007429-T1	0	Y	n2-9c17/18o
GPALN_005098-T1	0	Y	n4-11c16/17o
GPALN_004370-T1	0	Y	n3-13c18/19o
GPALN_015262-T1	0	Y	n31-44c49/50o
GPALN_012299-T1	0	Y	n11-21c26/27o
GPALN_014354-T1	0	Y	n10-21c25/26o
GPALN_007416-T1	0	Y	n10-21c26/27o
GPALN_009711-T1	0	Y	n3-14c19/20o
GPALN_003463-T1	0	Y	n2-13c18/19o
GPALN_001211-T1	0	Y	n2-7c11/12o
GPALN_000371-T1	0	Y	n4-15c27/28o
GPALN_014775-T1	0	Y	n5-15c20/21o
GPALN_007682-T1	0	Y	n3-13c18/19o
GPALN_015659-T1	0	Y	n4-16c20/21o
GPALN_016197-T1	0	Y	n4-15c20/21o
GPALN_006911-T1	0	Y	n4-15c23/24o
GPALN_008176-T1	0	Y	n3-14c19/20o
GPALN_007846-T1	1	Y	n3-18c24/25o325-341i
GPALN_002175-T1	0	Y	n3-14c19/20o
GPALN_003618-T1	0	Y	n3-10c15/16o
GPALN_004679-T1	0	Y	n4-14c19/20o
GPALN_014304-T2	0	Y	n8-23c28/29o
GPALN_000396-T1	0	Y	n13-21c26/27o
GPALN_016307-T1	0	Y	n4-14c22/23o
GPALN_003999-T1	0	Y	n5-12c20/21o
GPALN_016164-T1	0	Y	n4-19c28/29o
GPALN_015841-T1	0	Y	n11-19c31/32o
GPALN_015931-T2	11	Y	n7-18c27/28o43-62i74-92o98-122i134-155o175-195i271-292o320-342i363-384o390-412i424-445o465-491i
GPALN_002082-T1	0	Y	n6-17c24/25o
GPALN_004585-T1	0	Y	n15-27c32/33o
GPALN_008346-T1	1	Y	n2-9c15/16o43-69i
GPALN_007189-T1	1	Y	n8-16c22/23o170-188i
GPALN_011875-T1	0	Y	n8-19c23/24o
GPALN_014693-T1	1	Y	n4-15c23/24o47-72i
GPALN_007060-T1	0	Y	n8-20c24/25o
GPALN_003438-T1	0	Y	n7-18c23/24o
GPALN_008144-T1	0	Y	n2-13c20/21o
GPALN_003547-T1	0	Y	n15-26c34/35o
GPALN_010190-T1	0	Y	n5-18c23/24o
GPALN_010210-T1	0	Y	n3-10c18/19o
GPALN_009291-T1	0	Y	n12-22c27/28o
GPALN_002879-T1	0	Y	n2-12c17/18o
GPALN_009377-T1	1	Y	n10-21c39/40o549-571i
GPALN_009590-T1	0	Y	n8-19c23/24o
GPALN_004974-T1	0	Y	n6-17c21/22o
GPALN_004469-T1	1	Y	n6-18c23/24o506-528i
GPALN_014121-T1	0	Y	n2-12c17/18o
GPALN_007697-T1	0	Y	n3-13c18/19o
GPALN_008971-T1	2	Y	n3-10c18/19o150-169i190-209o
GPALN_002681-T1	0	Y	n11-23c28/29o
GPALN_007436-T1	0	Y	n4-14c24/25o
GPALN_013810-T1	0	Y	n7-18c26/27o
GPALN_008215-T1	0	Y	n9-16c21/22o
GPALN_000018-T1	5	Y	n7-18c22/23o527-545i565-590o610-635i944-961o967-986i
GPALN_010257-T1	0	Y	n3-14c19/20o
GPALN_014975-T1	0	Y	n2-12c20/21o
GPALN_002693-T1	0	Y	n14-22c27/28o
GPALN_004840-T2	0	Y	n7-17c22/23o
GPALN_004887-T1	0	Y	n11-23c29/30o
GPALN_015177-T1	0	Y	n8-19c24/25o
GPALN_007314-T1	0	Y	n3-14c18/19o
GPALN_002593-T1	0	Y	n8-18c22/23o
GPALN_012412-T1	0	Y	n3-11c16/17o
GPALN_005846-T1	0	Y	n3-16c24/25o
GPALN_008335-T2	0	Y	n8-19c23/24o
GPALN_014045-T1	0	Y	n5-16c21/22o
GPALN_007644-T3	0	Y	n13-23c28/29o
GPALN_003357-T1	0	Y	n6-19c24/25o
GPALN_010330-T1	0	Y	n12-23c28/29o
GPALN_004676-T1	1	Y	n7-22c27/28o215-237i
GPALN_000557-T1	0	Y	n10-21c25/26o
GPALN_012998-T2	0	Y	n6-16c23/24o
GPALN_000016-T1	0	Y	n2-12c17/18o
GPALN_004868-T1	0	Y	n7-17c22/23o
GPALN_004802-T1	0	Y	n2-13c19/20o
GPALN_013504-T2	0	Y	n2-9c13/14o
GPALN_004009-T1	0	Y	n4-19c24/25o
GPALN_007287-T1	0	Y	n3-21c29/30o
GPALN_012511-T1	0	Y	n8-17c21/22o
GPALN_002687-T1	0	Y	n14-25c29/30o
GPALN_004425-T1	0	Y	n2-13c17/18o
GPALN_015267-T1	0	Y	n8-18c23/24o
GPALN_013586-T1	0	Y	n4-15c20/21o
GPALN_001149-T1	0	Y	n4-15c24/25o
GPALN_005351-T1	0	Y	n3-11c15/16o
GPALN_004749-T1	0	Y	n5-16c23/24o
GPALN_007995-T1	0	Y	n4-17c21/22o
GPALN_012302-T1	1	Y	n7-17c24/25o270-292i
GPALN_008074-T1	0	Y	n3-18c23/24o
GPALN_016343-T1	0	Y	n4-17c22/23o
GPALN_006604-T1	0	Y	n6-17c25/26o
GPALN_015172-T1	0	Y	n3-14c19/20o
GPALN_001443-T1	0	Y	n3-11c16/17o
GPALN_001753-T1	0	Y	n3-10c14/15o
GPALN_004570-T1	0	Y	n15-25c33/34o
GPALN_010582-T1	0	Y	n7-19c25/26o
GPALN_002964-T1	0	Y	n8-16c22/23o
GPALN_011629-T1	0	Y	n5-13c18/19o
GPALN_003916-T1	0	Y	n27-38c42/43o
GPALN_010823-T2	2	Y	n18-29c34/35o50-72i408-427o
GPALN_004411-T1	0	Y	n5-15c23/24o
GPALN_003752-T1	0	Y	n3-11c19/20o
GPALN_007369-T1	0	Y	n6-16c24/25o
GPALN_015654-T2	0	Y	n10-21c26/27o
GPALN_014951-T1	0	Y	n8-16c21/22o
GPALN_003377-T1	0	Y	n3-16c24/25o
GPALN_004770-T1	0	Y	n3-15c19/20o
GPALN_005240-T1	1	Y	n5-15c21/22o52-75i
GPALN_011610-T1	0	Y	n2-12c17/18o
GPALN_010828-T1	0	Y	n4-19c28/29o
GPALN_012766-T1	0	Y	n6-16c22/23o
GPALN_010431-T1	0	Y	n7-16c21/22o
GPALN_001954-T1	0	Y	n6-17c22/23o
GPALN_011621-T1	0	Y	n5-13c18/19o
GPALN_005161-T1	0	Y	n3-14c19/20o
GPALN_010778-T1	0	Y	n6-16c22/23o
GPALN_002219-T1	0	Y	n3-14c26/27o
GPALN_004950-T1	0	Y	n4-12c22/23o
GPALN_012897-T1	0	Y	n3-13c18/19o
GPALN_013356-T1	0	Y	n3-14c21/22o
GPALN_013228-T1	1	Y	n2-12c16/17o145-167i
GPALN_009466-T1	0	Y	n8-16c23/24o
GPALN_000654-T1	0	Y	n4-16c21/22o
GPALN_014635-T1	0	Y	n13-24c32/33o
GPALN_007033-T1	0	Y	n4-15c24/25o
GPALN_001570-T1	0	Y	n2-10c18/19o
GPALN_005926-T1	0	Y	n2-20c25/26o
GPALN_007646-T1	0	Y	n13-24c29/30o
GPALN_013784-T1	1	Y	n16-27c32/33o251-277i
GPALN_011991-T1	1	Y	n10-18c23/24o33-55i
GPALN_013297-T1	0	Y	n8-15c20/21o
GPALN_014590-T1	0	Y	n2-12c20/21o
GPALN_014867-T1	0	Y	n5-16c21/22o
GPALN_015121-T1	0	Y	n6-17c24/25o
GPALN_007592-T1	0	Y	n3-15c20/21o
GPALN_009629-T1	0	Y	n8-19c24/25o
GPALN_015686-T1	0	Y	n7-14c19/20o
GPALN_015447-T1	0	Y	n14-25c30/31o
GPALN_005158-T1	1	Y	n4-12c17/18o644-673i
GPALN_010575-T1	0	Y	n7-19c24/25o
GPALN_005897-T1	1	Y	n7-18c22/23o332-350i
GPALN_009906-T1	0	Y	n7-17c22/23o
GPALN_000871-T1	1	Y	n2-12c19/20o29-49i
GPALN_007708-T1	2	Y	n3-18c23/24o637-655i662-685o
GPALN_011396-T1	0	Y	n27-37c46/47o
GPALN_007295-T1	0	Y	n3-21c26/27o
GPALN_004854-T1	0	Y	n4-16c21/22o
GPALN_011593-T1	0	Y	n22-30c48/49o
GPALN_013749-T1	0	Y	n4-15c26/27o
GPALN_003264-T1	0	Y	n3-13c17/18o
GPALN_014960-T1	0	Y	n6-17c25/26o
GPALN_002490-T1	1	Y	n7-19c24/25o838-858i
GPALN_016000-T1	1	Y	n6-17c22/23o226-251i
GPALN_000948-T1	0	Y	n12-23c28/29o
GPALN_008859-T1	0	Y	n2-13c18/19o
GPALN_009907-T1	0	Y	n7-17c22/23o
GPALN_002314-T1	0	Y	n8-20c24/25o
GPALN_004546-T1	0	Y	n18-29c38/39o
GPALN_003828-T1	0	Y	n4-15c20/21o
GPALN_005765-T1	0	Y	n7-17c22/23o
GPALN_016089-T1	0	Y	n6-15c20/21o
GPALN_006587-T1	0	Y	n9-20c25/26o
GPALN_009566-T1	1	Y	n4-11c16/17o334-354i
GPALN_007380-T1	0	Y	n4-15c27/28o
GPALN_002472-T1	9	Y	n16-24c31/32o796-816i881-901o907-930i951-969o975-991i1012-1031o1078-1097i1109-1130o1136-1157i
GPALN_013423-T1	0	Y	n4-15c20/21o
GPALN_015285-T1	0	Y	n5-16c21/22o
GPALN_006406-T1	4	Y	n10-21c30/31o192-218i281-306o318-344i1304-1326o
GPALN_016360-T1	0	Y	n10-20c24/25o
GPALN_009608-T1	0	Y	n15-26c30/31o
GPALN_013168-T1	0	Y	n5-20c24/25o
GPALN_015701-T1	0	Y	n5-16c20/21o
GPALN_004830-T1	0	Y	n3-14c19/20o
GPALN_014546-T1	0	Y	n10-20c25/26o
GPALN_012803-T2	0	Y	n18-29c33/34o
GPALN_009512-T1	1	Y	n3-15c19/20o355-373i
GPALN_003426-T1	0	Y	n9-17c22/23o
GPALN_007672-T1	0	Y	n6-19c27/28o
GPALN_008277-T1	1	Y	n5-14c19/20o1020-1041i
GPALN_015921-T1	0	Y	n2-10c15/16o
GPALN_004712-T1	0	Y	n7-18c27/28o
GPALN_005896-T1	0	Y	n3-10c15/16o
GPALN_002707-T1	0	Y	n7-18c26/27o
GPALN_011197-T1	6	Y	n3-14c19/20o1280-1302i1348-1368o1487-1513i1534-1558o1573-1594i1826-1847o
GPALN_009441-T1	0	Y	n5-18c23/24o
GPALN_010973-T1	1	Y	n9-19c26/27o180-201i
GPALN_002769-T1	0	Y	n7-18c26/27o
GPALN_013481-T1	0	Y	n3-14c22/23o
GPALN_008444-T1	0	Y	n19-30c35/36o
GPALN_014107-T1	0	Y	n7-20c25/26o
GPALN_014086-T1	0	Y	n4-13c21/22o
GPALN_011614-T1	0	Y	n3-13c18/19o
GPALN_014594-T1	0	Y	n3-10c15/16o
GPALN_004479-T1	1	Y	n6-18c22/23o68-86i
GPALN_004592-T1	0	Y	n131-142c154/155o
GPALN_012912-T1	0	Y	n10-21c27/28o
GPALN_010424-T1	0	Y	n3-13c17/18o
GPALN_002969-T1	0	Y	n5-20c24/25o
GPALN_005120-T1	0	Y	n10-20c28/29o
GPALN_005745-T1	1	Y	n3-14c19/20o289-308i
GPALN_016289-T1	0	Y	n7-15c19/20o
GPALN_005349-T1	0	Y	n4-14c23/24o
GPALN_010904-T1	0	Y	n5-16c28/29o
GPALN_001845-T2	0	Y	n3-14c19/20o
GPALN_006215-T2	0	Y	n4-13c18/19o
GPALN_007194-T1	1	Y	n3-14c18/19o226-243i
GPALN_013421-T1	0	Y	n9-19c24/25o
GPALN_012478-T1	0	Y	n3-13c18/19o
GPALN_010570-T1	0	Y	n3-11c16/17o
GPALN_004759-T1	1	Y	n9-17c25/26o35-53i
GPALN_004179-T1	1	Y	n3-14c19/20o291-311i
GPALN_009147-T1	0	Y	n8-19c23/24o
GPALN_010171-T1	0	Y	n3-14c19/20o
GPALN_009641-T1	0	Y	n6-17c22/23o
GPALN_003430-T1	0	Y	n7-17c22/23o
GPALN_008814-T1	1	Y	n17-28c33/34o110-128i
GPALN_010842-T1	1	Y	n15-30c35/36o51-72i
GPALN_010599-T1	0	Y	n5-18c26/27o
GPALN_006650-T1	0	Y	n3-15c20/21o
GPALN_015408-T1	0	Y	n6-16c22/23o
GPALN_003953-T1	0	Y	n10-21c26/27o
GPALN_013222-T1	0	Y	n4-11c16/17o
GPALN_015227-T1	0	Y	n2-10c17/18o
GPALN_009716-T1	1	Y	n10-21c28/29o154-177i
GPALN_014881-T1	0	Y	n4-15c24/25o
GPALN_007357-T1	0	Y	n5-15c24/25o
GPALN_002419-T1	1	Y	n14-24c31/32o237-260i
GPALN_000771-T1	0	Y	n19-29c37/38o
GPALN_002163-T1	0	Y	n2-11c16/17o
GPALN_013277-T1	0	Y	n6-16c21/22o
GPALN_008167-T1	0	Y	n2-10c16/17o
GPALN_003434-T1	0	Y	n16-27c32/33o
GPALN_003606-T1	0	Y	n10-21c25/26o
GPALN_014379-T1	0	Y	n3-14c19/20o
GPALN_008417-T1	0	Y	n17-25c31/32o
GPALN_007467-T1	0	Y	n11-23c27/28o
GPALN_002890-T1	0	Y	n27-38c43/44o
GPALN_000414-T1	0	Y	n8-20c25/26o
GPALN_008680-T1	0	Y	n8-22c26/27o
GPALN_007445-T1	0	Y	n4-19c23/24o
GPALN_002334-T1	0	Y	n3-13c17/18o
GPALN_007886-T1	1	Y	n12-23c29/30o197-222i
GPALN_007385-T1	0	Y	n4-16c25/26o
GPALN_001866-T1	0	Y	n4-12c17/18o
GPALN_015622-T2	0	Y	n2-13c20/21o
GPALN_014575-T1	0	Y	n2-13c17/18o
GPALN_002641-T1	6	Y	n6-17c22/23o587-612i624-644o697-718i739-763o775-798i810-832o
GPALN_007541-T1	0	Y	n4-15c23/24o
GPALN_004117-T1	4	Y	n10-23c28/29o185-206i248-268o274-292i304-324o
GPALN_004131-T1	0	Y	n2-13c17/18o
GPALN_008233-T1	0	Y	n10-21c26/27o
GPALN_000433-T1	0	Y	n7-18c23/24o
GPALN_001617-T1	0	Y	n3-14c19/20o
GPALN_012894-T1	0	Y	n4-14c19/20o
GPALN_011587-T2	0	Y	n12-23c28/29o
GPALN_015757-T1	4	Y	n16-23c27/28o200-223i254-275o281-299i311-331o
GPALN_010285-T1	0	Y	n8-16c21/22o
GPALN_010403-T1	0	Y	n3-13c17/18o
GPALN_000936-T1	0	Y	n6-17c25/26o
GPALN_008619-T1	0	Y	n4-15c23/24o
GPALN_002150-T1	0	Y	n3-14c19/20o
GPALN_013747-T1	1	Y	n14-27c32/33o484-507i
GPALN_014381-T1	0	Y	n3-14c19/20o
GPALN_011521-T1	0	Y	n8-16c21/22o
GPALN_012004-T1	0	Y	n7-18c26/27o
GPALN_006286-T1	0	Y	n2-12c17/18o
GPALN_010180-T1	0	Y	n8-17c21/22o
GPALN_015137-T1	0	Y	n2-9c13/14o
GPALN_003010-T1	0	Y	n6-18c23/24o
GPALN_013082-T1	0	Y	n2-13c18/19o
GPALN_002365-T1	0	Y	n5-15c20/21o
GPALN_009443-T1	0	Y	n9-19c27/28o
GPALN_013824-T1	1	Y	n17-30c34/35o207-230i
GPALN_009492-T1	0	Y	n8-18c22/23o
GPALN_010521-T1	0	Y	n7-18c25/26o
GPALN_004729-T1	0	Y	n6-17c22/23o
GPALN_008420-T1	0	Y	n2-16c21/22o
GPALN_015428-T1	1	Y	n19-30c34/35o252-270i
GPALN_014324-T1	0	Y	n4-15c23/24o
GPALN_009825-T1	0	Y	n5-15c27/28o
GPALN_004747-T1	0	Y	n6-13c18/19o
GPALN_007130-T1	0	Y	n6-16c21/22o
GPALN_010199-T1	0	Y	n2-12c17/18o
GPALN_009626-T1	0	Y	n33-44c49/50o
GPALN_002329-T1	0	Y	n2-10c14/15o
GPALN_001317-T1	0	Y	n3-15c20/21o
GPALN_007372-T1	0	Y	n6-18c28/29o
GPALN_004014-T1	0	Y	n7-18c22/23o
GPALN_013818-T1	12	Y	n36-49c53/54o567-586i598-620o672-690i702-721o727-750i959-975o1015-1036i1048-1065o1071-1094i1115-1137o1152-1172i1184-1201o
GPALN_010513-T1	0	Y	n4-16c21/22o
GPALN_016277-T1	1	Y	n3-14c19/20o241-261i
GPALN_004092-T1	0	Y	n10-21c33/34o
GPALN_009638-T1	0	Y	n6-17c22/23o
GPALN_007549-T1	0	Y	n3-13c20/21o
GPALN_004376-T2	0	Y	n8-23c31/32o
GPALN_003306-T1	0	Y	n3-16c21/22o
GPALN_002819-T1	0	Y	n36-48c53/54o
GPALN_010313-T1	0	Y	n6-14c22/23o
GPALN_009786-T1	0	Y	n3-13c17/18o
GPALN_015243-T1	0	Y	n8-19c24/25o
GPALN_004856-T1	0	Y	n9-20c25/26o
GPALN_004018-T1	0	Y	n6-16c25/26o
GPALN_015799-T1	0	Y	n3-14c18/19o
GPALN_010801-T1	0	Y	n4-11c19/20o
GPALN_012624-T1	2	Y	n3-14c22/23o858-879i1138-1156o
GPALN_002018-T2	0	Y	n2-13c25/26o
GPALN_015418-T1	0	Y	n9-20c28/29o
GPALN_014880-T1	0	Y	n5-17c22/23o
GPALN_013990-T1	0	Y	n7-18c23/24o
GPALN_005905-T1	0	Y	n6-17c22/23o
GPALN_011858-T1	0	Y	n4-15c20/21o
GPALN_012818-T1	0	Y	n5-20c28/29o
GPALN_004019-T1	0	Y	n3-14c19/20o
GPALN_005137-T1	1	Y	n4-15c23/24o263-286i
GPALN_006730-T1	0	Y	n9-16c21/22o
GPALN_005082-T1	0	Y	n8-19c24/25o
GPALN_007502-T1	0	Y	n4-14c18/19o
GPALN_008339-T1	1	Y	n7-20c25/26o146-165i
GPALN_015334-T1	0	Y	n7-16c20/21o
GPALN_001318-T1	0	Y	n7-15c20/21o
GPALN_006754-T1	0	Y	n16-29c34/35o
GPALN_002098-T1	0	Y	n7-18c23/24o
GPALN_004897-T1	0	Y	n4-15c20/21o
GPALN_015287-T1	1	Y	n4-14c19/20o247-268i
GPALN_016341-T1	0	Y	n8-19c26/27o
GPALN_011230-T1	0	Y	n4-11c16/17o
GPALN_001944-T1	0	Y	n11-22c26/27o
GPALN_002364-T1	0	Y	n4-13c21/22o
GPALN_014862-T1	0	Y	n6-16c20/21o
GPALN_013081-T1	0	Y	n5-14c24/25o
GPALN_001449-T1	0	Y	n2-11c16/17o
GPALN_010360-T1	0	Y	n9-20c25/26o
GPALN_007737-T1	1	Y	n8-19c24/25o112-136i
GPALN_004294-T1	0	Y	n6-17c25/26o
GPALN_002574-T1	0	Y	n5-13c18/19o
GPALN_006041-T1	2	Y	n10-21c26/27o75-95i107-124o
GPALN_011956-T1	0	Y	n2-10c15/16o
GPALN_009103-T1	1	Y	n11-22c31/32o312-334i
GPALN_002775-T1	0	Y	n8-19c25/26o
GPALN_003770-T1	1	Y	n7-15c20/21o148-170i
GPALN_008646-T1	0	Y	n8-15c20/21o
GPALN_004174-T1	0	Y	n3-14c26/27o
GPALN_005489-T1	0	Y	n4-15c20/21o
GPALN_007837-T1	0	Y	n3-13c24/25o
GPALN_003383-T1	0	Y	n5-15c20/21o
GPALN_012915-T1	0	Y	n4-15c33/34o
GPALN_006126-T1	0	Y	n9-19c26/27o
GPALN_006719-T1	0	Y	n7-18c24/25o
GPALN_012298-T1	0	Y	n11-21c26/27o
GPALN_001478-T1	0	Y	n5-15c20/21o
GPALN_003786-T1	0	Y	n8-19c24/25o
GPALN_011532-T1	0	Y	n4-15c20/21o
GPALN_001216-T1	0	Y	n4-15c20/21o
GPALN_010586-T1	0	Y	n15-25c29/30o
GPALN_015163-T1	0	Y	n2-10c15/16o
GPALN_003693-T1	0	Y	n7-18c24/25o
GPALN_005084-T1	0	Y	n2-10c15/16o
GPALN_002068-T1	0	Y	n6-13c17/18o
GPALN_011097-T1	0	Y	n7-18c23/24o
GPALN_014800-T1	0	Y	n6-16c21/22o
GPALN_014235-T1	0	Y	n4-15c23/24o
GPALN_004771-T1	11	Y	n18-29c34/35o44-66i78-98o127-150i162-186o192-212i268-292o304-327i334-352o358-382i394-418o478-498i
GPALN_004597-T1	0	Y	n4-16c21/22o
GPALN_002204-T1	0	Y	n8-19c26/27o
GPALN_011319-T2	0	Y	n5-12c17/18o
GPALN_010019-T1	0	Y	n4-11c16/17o
GPALN_009839-T1	0	Y	n6-17c22/23o
GPALN_005776-T1	0	Y	n2-13c18/19o
GPALN_015742-T1	1	Y	n4-15c19/20o716-734i
GPALN_007805-T1	0	Y	n2-13c20/21o
GPALN_005074-T1	0	Y	n157-162c166/167o
GPALN_002710-T1	0	Y	n9-20c25/26o
GPALN_001950-T2	5	Y	n4-12c17/18o33-52i64-81o93-113i122-142o162-181i
GPALN_015230-T1	0	Y	n3-15c23/24o
GPALN_012732-T1	0	Y	n15-23c28/29o
GPALN_000198-T1	0	Y	n19-30c35/36o
GPALN_007947-T1	0	Y	n8-16c21/22o
GPALN_009357-T1	0	Y	n4-14c19/20o
GPALN_008104-T1	0	Y	n2-10c15/16o
GPALN_001668-T1	0	Y	n10-21c29/30o
GPALN_005032-T1	0	Y	n2-10c18/19o
GPALN_013609-T1	2	Y	n4-17c22/23o32-52i239-256o
GPALN_010724-T1	0	Y	n6-17c22/23o
GPALN_002387-T1	0	Y	n6-17c22/23o
GPALN_012007-T1	0	Y	n4-19c23/24o
GPALN_016083-T1	0	Y	n4-16c21/22o
GPALN_013923-T1	0	Y	n5-18c23/24o
GPALN_010883-T1	0	Y	n11-21c29/30o
GPALN_013261-T1	0	Y	n3-13c18/19o
GPALN_015578-T1	0	Y	n9-20c25/26o
GPALN_010937-T1	0	Y	n16-23c28/29o
GPALN_015238-T1	0	Y	n3-13c21/22o
GPALN_016330-T1	0	Y	n4-18c24/25o
GPALN_016199-T1	0	Y	n5-17c25/26o
GPALN_011715-T1	0	Y	n5-16c22/23o
GPALN_002399-T1	0	Y	n15-24c29/30o
GPALN_010600-T1	0	Y	n5-16c24/25o
GPALN_003431-T1	0	Y	n8-18c23/24o
GPALN_002920-T1	0	Y	n2-13c18/19o
GPALN_005145-T1	0	Y	n3-15c19/20o
GPALN_005809-T1	0	Y	n4-15c24/25o
GPALN_014904-T1	0	Y	n14-25c32/33o
GPALN_016166-T1	0	Y	n6-19c28/29o
GPALN_015545-T1	0	Y	n6-14c19/20o
GPALN_013577-T1	0	Y	n8-20c26/27o
GPALN_012589-T1	0	Y	n8-19c24/25o
GPALN_007696-T1	0	Y	n3-14c19/20o
GPALN_002893-T1	0	Y	n3-14c19/20o
GPALN_001548-T1	0	Y	n2-10c17/18o
GPALN_011522-T1	0	Y	n8-16c23/24o
GPALN_000933-T1	3	Y	n16-26c30/31o40-62i83-101o113-135i
GPALN_010823-T1	4	Y	n10-19c24/25o252-274i281-298o318-339i675-694o
GPALN_002383-T1	0	Y	n5-15c23/24o
GPALN_009355-T1	2	Y	n16-27c32/33o69-88i100-120o
GPALN_007823-T1	0	Y	n14-21c26/27o
GPALN_010160-T1	1	Y	n3-11c16/17o306-324i
GPALN_013276-T1	0	Y	n2-12c17/18o
GPALN_007550-T1	0	Y	n7-15c20/21o
GPALN_009637-T1	0	Y	n6-17c22/23o
GPALN_008620-T1	0	Y	n3-14c26/27o
GPALN_015822-T1	0	Y	n3-11c16/17o
GPALN_012634-T1	0	Y	n34-45c52/53o
GPALN_005167-T1	0	Y	n3-15c23/24o
GPALN_011627-T1	0	Y	n4-17c25/26o
GPALN_012455-T1	1	Y	n8-19c26/27o36-54i
GPALN_011464-T1	0	Y	n6-17c25/26o
GPALN_003090-T1	0	Y	n8-20c25/26o
GPALN_000560-T1	3	Y	n4-11c15/16o25-43i55-73o130-152i
GPALN_015359-T1	1	Y	n7-19c25/26o351-372i
GPALN_001991-T1	0	Y	n5-13c18/19o
GPALN_012022-T1	0	Y	n8-20c28/29o
GPALN_012532-T1	0	Y	n13-23c28/29o
GPALN_003820-T1	0	Y	n4-15c23/24o
GPALN_004076-T1	0	Y	n2-13c17/18o
GPALN_013004-T1	0	Y	n5-15c23/24o
GPALN_010720-T1	0	Y	n7-17c22/23o
GPALN_005418-T1	0	Y	n15-25c31/32o
GPALN_009908-T1	0	Y	n7-17c22/23o
GPALN_001773-T1	0	Y	n6-17c22/23o
GPALN_010679-T1	11	Y	n15-26c36/37o125-145i152-171o177-199i211-234o246-266i314-336o356-377i389-410o416-440i447-469o481-505i
GPALN_010628-T2	0	Y	n14-24c32/33o
GPALN_015155-T1	0	Y	n4-15c22/23o
GPALN_004667-T1	0	Y	n20-30c35/36o
GPALN_006588-T1	0	Y	n2-13c18/19o
GPALN_005751-T1	0	Y	n3-14c20/21o
GPALN_004763-T1	0	Y	n7-17c25/26o
GPALN_004536-T1	9	Y	n12-23c28/29o209-229i296-315o321-342i363-382o388-404i425-444o490-507i519-539o551-572i
GPALN_009586-T1	0	Y	n5-16c24/25o
GPALN_007176-T1	1	Y	n2-13c18/19o239-260i
GPALN_008289-T1	1	Y	n14-25c32/33o237-262i
GPALN_012925-T1	0	Y	n16-27c35/36o
GPALN_012702-T1	1	Y	n3-14c26/27o370-388i
GPALN_006452-T1	0	Y	n7-15c24/25o
GPALN_014032-T1	0	Y	n9-20c24/25o
GPALN_008558-T1	0	Y	n6-17c21/22o
GPALN_003486-T1	0	Y	n3-14c19/20o
GPALN_004476-T1	3	Y	n2-13c17/18o27-44i56-74o80-98i
GPALN_003983-T1	1	Y	n2-13c18/19o394-412i
GPALN_001070-T1	1	Y	n8-19c24/25o434-454i
GPALN_015296-T1	0	Y	n4-15c20/21o
GPALN_010591-T1	0	Y	n4-14c26/27o
GPALN_010509-T1	0	Y	n3-14c19/20o
GPALN_005113-T1	0	Y	n6-17c29/30o
GPALN_006802-T1	0	Y	n6-16c23/24o
GPALN_010295-T1	0	Y	n6-17c22/23o
GPALN_010643-T1	0	Y	n3-14c19/20o
GPALN_012515-T1	0	Y	n4-16c21/22o
GPALN_006769-T1	0	Y	n16-29c34/35o
GPALN_000383-T1	0	Y	n2-17c22/23o
GPALN_011838-T1	0	Y	n3-13c18/19o
GPALN_008779-T1	0	Y	n21-31c38/39o
GPALN_007884-T1	0	Y	n12-22c30/31o
GPALN_006977-T1	0	Y	n2-13c18/19o
GPALN_001004-T1	0	Y	n4-17c22/23o
GPALN_014639-T1	0	Y	n4-23c28/29o
GPALN_012647-T1	0	Y	n2-12c17/18o
GPALN_015233-T1	0	Y	n2-12c17/18o
GPALN_010625-T1	0	Y	n13-25c30/31o
GPALN_010148-T1	0	Y	n6-16c20/21o
GPALN_013246-T1	0	Y	n9-19c27/28o
GPALN_015286-T1	1	Y	n4-14c19/20o250-271i
GPALN_006430-T1	0	Y	n25-36c44/45o
GPALN_012664-T1	1	Y	n11-26c31/32o577-593i
GPALN_005769-T1	0	Y	n4-15c20/21o
GPALN_008783-T1	0	Y	n8-23c27/28o
GPALN_010603-T1	0	Y	n2-13c18/19o
GPALN_004203-T1	0	Y	n2-13c20/21o
GPALN_006407-T1	1	Y	n16-27c33/34o197-221i
GPALN_005180-T2	0	Y	n3-14c28/29o
GPALN_008308-T1	1	Y	n2-10c15/16o753-772i
GPALN_010232-T1	0	Y	n4-15c20/21o
GPALN_004798-T1	0	Y	n7-17c25/26o
GPALN_000510-T1	0	Y	n4-15c20/21o
GPALN_003803-T1	2	Y	n9-21c26/27o78-98i105-122o
GPALN_014479-T1	0	Y	n6-13c18/19o
GPALN_003340-T1	0	Y	n6-17c25/26o
GPALN_012344-T1	1	Y	n34-43c48/49o58-79i
GPALN_000705-T1	0	Y	n6-16c23/24o
GPALN_002861-T1	0	Y	n2-9c14/15o
GPALN_008377-T1	0	Y	n4-15c22/23o
GPALN_011705-T1	0	Y	n14-25c32/33o
GPALN_005209-T2	2	Y	n10-22c27/28o156-177i189-209o
GPALN_000452-T1	0	Y	n17-28c33/34o
GPALN_014102-T1	0	Y	n8-19c24/25o
GPALN_001165-T1	0	Y	n3-11c16/17o
GPALN_001525-T1	1	Y	n2-17c21/22o233-249i
GPALN_011558-T1	1	Y	n7-15c23/24o342-360i
GPALN_014327-T1	0	Y	n4-15c23/24o
GPALN_001640-T1	1	Y	n17-28c35/36o175-198i
GPALN_009898-T1	0	Y	n12-23c30/31o
GPALN_001161-T1	0	Y	n8-19c23/24o
GPALN_006035-T1	2	Y	n10-21c25/26o572-589i1000-1018o
GPALN_002682-T1	5	Y	n8-19c24/25o288-310i317-338o391-409i421-446o500-519i
GPALN_015804-T1	0	Y	n3-14c18/19o
GPALN_015733-T1	0	Y	n8-26c32/33o
GPALN_012703-T1	0	Y	n7-18c24/25o
GPALN_002802-T1	0	Y	n7-22c28/29o
GPALN_016023-T1	1	Y	n4-15c20/21o60-85i
GPALN_001147-T1	1	Y	n6-21c26/27o131-149i
GPALN_012122-T1	1	Y	n7-25c30/31o244-262i
GPALN_007811-T1	0	Y	n6-17c25/26o
GPALN_004007-T1	0	Y	n3-14c19/20o
GPALN_014103-T1	0	Y	n52-59c67/68o
GPALN_012374-T1	1	Y	n4-11c16/17o40-60i
GPALN_002844-T1	1	Y	n3-13c18/19o271-292i
GPALN_015629-T1	0	Y	n2-13c20/21o
GPALN_008449-T1	0	Y	n28-39c44/45o
GPALN_008500-T1	0	Y	n4-14c19/20o
GPALN_001922-T1	9	Y	n5-16c20/21o30-45i50-68o74-91i103-127o147-168i257-274o280-297i330-349o369-392i
GPALN_016244-T1	0	Y	n6-13c18/19o
GPALN_003375-T1	0	Y	n6-19c31/32o
GPALN_002546-T1	0	Y	n6-18c23/24o
GPALN_016368-T2	0	Y	n7-15c19/20o
GPALN_009634-T1	1	Y	n16-27c32/33o280-298i
GPALN_001951-T1	0	Y	n4-18c23/24o
GPALN_001150-T1	0	Y	n4-15c23/24o
GPALN_014034-T1	0	Y	n3-14c26/27o
GPALN_011455-T1	0	Y	n5-12c17/18o
GPALN_004172-T1	0	Y	n3-15c20/21o
GPALN_016117-T1	0	Y	n6-17c22/23o
GPALN_001151-T1	0	Y	n6-21c26/27o
GPALN_003797-T1	0	Y	n8-17c21/22o
GPALN_010561-T1	0	Y	n26-37c42/43o
GPALN_006586-T1	0	Y	n7-18c26/27o
GPALN_012815-T1	0	Y	n9-18c23/24o
GPALN_008232-T1	0	Y	n5-16c25/26o
GPALN_007544-T1	0	Y	n5-15c19/20o
GPALN_012152-T1	0	Y	n7-15c20/21o
GPALN_006088-T1	1	Y	n3-15c23/24o237-254i
GPALN_000866-T1	1	Y	n12-25c30/31o54-73i
GPALN_011376-T1	0	Y	n7-22c27/28o
GPALN_002442-T1	0	Y	n17-27c32/33o
GPALN_000994-T1	0	Y	n5-16c23/24o
GPALN_003809-T1	0	Y	n4-15c20/21o
GPALN_009860-T1	0	Y	n10-21c26/27o
GPALN_009431-T1	0	Y	n3-14c19/20o
GPALN_007464-T1	0	Y	n12-19c27/28o
GPALN_007229-T1	0	Y	n9-20c38/39o
GPALN_002542-T1	1	Y	n3-14c18/19o446-467i
GPALN_004293-T1	0	Y	n5-16c34/35o
GPALN_005870-T1	0	Y	n6-17c22/23o
GPALN_014876-T1	0	Y	n2-13c17/18o
GPALN_013613-T1	1	Y	n14-21c29/30o612-635i
GPALN_016098-T2	0	Y	n6-17c24/25o
GPALN_007638-T1	0	Y	n3-14c19/20o
GPALN_016098-T1	1	Y	n6-17c24/25o135-161i
GPALN_001091-T1	0	Y	n3-14c19/20o
GPALN_010449-T1	0	Y	n15-26c31/32o
GPALN_010809-T1	1	Y	n3-14c18/19o477-497i
GPALN_010319-T1	0	Y	n7-18c23/24o
GPALN_013399-T1	0	Y	n5-16c22/23o
GPALN_002252-T2	0	Y	n4-15c21/22o
GPALN_009111-T1	0	Y	n11-22c34/35o
GPALN_006759-T1	0	Y	n11-22c27/28o
GPALN_003719-T1	0	Y	n7-17c22/23o
GPALN_005882-T1	0	Y	n6-14c19/20o
GPALN_015209-T1	0	Y	n9-16c23/24o
GPALN_007070-T1	0	Y	n4-15c20/21o
GPALN_015873-T1	4	Y	n6-17c22/23o46-69i89-109o142-160i350-371o
GPALN_003334-T1	0	Y	n8-21c26/27o
GPALN_009406-T1	0	Y	n4-22c28/29o
GPALN_014077-T2	0	Y	n15-26c34/35o
GPALN_002509-T1	0	Y	n5-16c24/25o
GPALN_002947-T1	0	Y	n4-18c22/23o
GPALN_014033-T1	0	Y	n3-12c16/17o
GPALN_003557-T1	1	Y	n7-19c23/24o194-211i
GPALN_014350-T1	0	Y	n5-16c24/25o
GPALN_008803-T1	6	Y	n17-27c31/32o923-940i1200-1228o1248-1272i1284-1308o1320-1341i1362-1384o
GPALN_014879-T1	0	Y	n2-13c18/19o
GPALN_010663-T1	0	Y	n18-29c34/35o
GPALN_013239-T1	0	Y	n2-12c16/17o
GPALN_006490-T1	3	Y	n2-10c15/16o334-357i369-390o402-424i
GPALN_008083-T1	0	Y	n9-22c29/30o
GPALN_014549-T1	0	Y	n3-14c18/19o
GPALN_000179-T2	0	Y	n15-25c29/30o
GPALN_005230-T1	1	Y	n12-23c29/30o198-217i
GPALN_006128-T1	1	Y	n11-21c26/27o50-68i
GPALN_012067-T1	0	Y	n4-19c24/25o
GPALN_004450-T1	0	Y	n3-13c22/23o
GPALN_003177-T1	0	Y	n4-15c23/24o
GPALN_003402-T1	0	Y	n4-16c21/22o
GPALN_016189-T1	0	Y	n3-14c19/20o
GPALN_001315-T1	0	Y	n5-13c19/20o
GPALN_004445-T1	0	Y	n3-13c22/23o
GPALN_005679-T1	9	Y	n19-30c35/36o292-313i325-346o506-527i548-574o1004-1025i1084-1102o1122-1138i1150-1174o1186-1209i
GPALN_010912-T1	1	Y	n5-16c20/21o68-88i
GPALN_008673-T1	1	Y	n19-30c42/43o1028-1052i
GPALN_015346-T1	0	Y	n12-25c31/32o
GPALN_007219-T1	0	Y	n7-18c23/24o
GPALN_000125-T1	0	Y	n4-12c17/18o
GPALN_003040-T1	0	Y	n3-14c19/20o
GPALN_000153-T1	0	Y	n20-33c38/39o
GPALN_005724-T1	1	Y	n7-15c20/21o142-165i
GPALN_002478-T1	0	Y	n7-18c26/27o
GPALN_006072-T1	0	Y	n6-20c25/26o
GPALN_008088-T1	0	Y	n12-23c31/32o
GPALN_009000-T1	0	Y	n21-31c36/37o
GPALN_010128-T1	0	Y	n4-15c20/21o
GPALN_007513-T1	0	Y	n4-12c21/22o
GPALN_010510-T1	0	Y	n6-17c23/24o
GPALN_004709-T1	0	Y	n21-32c37/38o
GPALN_010859-T1	0	Y	n53-59c63/64o
GPALN_012983-T1	4	Y	n4-15c20/21o30-47i127-144o150-170i238-264o
GPALN_007375-T1	0	Y	n3-11c19/20o
GPALN_002432-T1	0	Y	n18-27c31/32o
GPALN_012124-T1	0	Y	n3-11c16/17o
GPALN_011615-T1	0	Y	n4-14c18/19o
GPALN_007976-T1	0	Y	n6-16c26/27o
GPALN_011972-T1	0	Y	n5-15c19/20o
GPALN_009823-T1	0	Y	n7-18c23/24o
GPALN_011472-T1	0	Y	n6-14c19/20o
GPALN_011921-T1	0	Y	n6-16c20/21o
GPALN_008303-T1	0	Y	n5-16c20/21o
GPALN_004498-T1	0	Y	n30-42c46/47o
GPALN_007919-T1	0	Y	n5-16c20/21o
GPALN_014094-T1	0	Y	n4-19c24/25o
GPALN_010144-T2	5	Y	n7-14c19/20o56-75i96-119o139-166i187-211o223-241i
GPALN_007193-T1	0	Y	n4-13c18/19o
GPALN_002382-T1	0	Y	n9-20c25/26o
GPALN_015223-T1	0	Y	n9-20c32/33o
GPALN_007186-T1	0	Y	n4-15c20/21o
GPALN_000742-T1	0	Y	n6-25c30/31o
GPALN_007269-T1	0	Y	n2-12c17/18o
GPALN_012873-T1	0	Y	n19-30c37/38o
GPALN_010632-T1	0	Y	n9-18c23/24o
GPALN_016050-T1	1	Y	n11-23c28/29o227-247i
GPALN_005452-T1	0	Y	n8-19c31/32o
GPALN_005743-T1	2	Y	n2-12c16/17o191-207i228-248o
GPALN_014967-T1	0	Y	n5-16c25/26o
GPALN_005117-T1	1	Y	n11-21c26/27o906-925i
GPALN_013682-T1	0	Y	n3-14c19/20o
GPALN_003905-T1	0	Y	n4-14c19/20o
GPALN_001369-T1	0	Y	n10-21c25/26o
GPALN_013064-T1	0	Y	n15-26c34/35o
GPALN_009837-T1	0	Y	n19-30c35/36o
GPALN_005273-T1	12	Y	n5-16c24/25o157-176i358-375o427-446i453-472o478-498i550-575o604-628i836-855o861-882i889-910o930-951i963-986o
GPALN_010588-T1	0	Y	n5-16c21/22o
GPALN_015324-T1	1	Y	n8-19c27/28o167-190i
GPALN_004275-T1	0	Y	n11-19c27/28o
GPALN_006038-T1	0	Y	n5-16c21/22o
GPALN_009905-T1	0	Y	n7-14c19/20o
GPALN_014334-T1	0	Y	n5-16c24/25o
GPALN_002252-T3	0	Y	n4-15c21/22o
GPALN_006059-T1	0	Y	n3-14c23/24o
GPALN_010067-T1	0	Y	n12-23c31/32o
GPALN_006026-T1	0	Y	n3-12c17/18o
GPALN_010290-T1	0	Y	n5-16c21/22o
GPALN_013594-T1	0	Y	n4-17c21/22o
GPALN_012416-T1	0	Y	n5-15c20/21o
GPALN_003757-T2	0	Y	n3-14c19/20o
GPALN_006022-T1	0	Y	n3-12c17/18o
GPALN_007899-T1	1	Y	n12-22c30/31o217-244i
GPALN_010614-T1	0	Y	n5-16c21/22o
GPALN_008594-T3	0	Y	n9-14c19/20o
GPALN_008745-T1	1	Y	n7-19c24/25o380-400i
GPALN_007233-T1	0	Y	n15-26c32/33o
GPALN_006424-T1	0	Y	n20-33c37/38o
GPALN_003588-T1	0	Y	n5-18c22/23o
GPALN_003077-T1	0	Y	n10-20c25/26o
GPALN_004065-T1	0	Y	n7-18c23/24o
GPALN_007202-T1	1	Y	n5-15c19/20o193-211i
GPALN_002988-T1	0	Y	n6-14c18/19o
GPALN_001591-T1	1	Y	n9-16c20/21o108-129i
GPALN_010542-T1	0	Y	n4-14c21/22o
GPALN_002987-T1	1	Y	n8-16c22/23o240-261i
GPALN_008379-T1	0	Y	n4-19c23/24o
GPALN_010571-T1	0	Y	n3-16c24/25o
GPALN_004678-T1	0	Y	n2-11c16/17o
GPALN_001950-T1	9	Y	n4-12c17/18o33-52i64-81o93-113i122-142o162-182i203-221o241-262i269-290o296-316i
GPALN_007315-T1	0	Y	n7-18c23/24o
GPALN_013273-T1	0	Y	n4-14c26/27o
GPALN_012409-T1	0	Y	n3-14c20/21o
GPALN_005712-T1	1	Y	n5-16c24/25o365-382i
GPALN_006176-T1	0	Y	n21-32c40/41o
GPALN_009886-T1	0	Y	n2-10c17/18o
GPALN_002194-T1	0	Y	n6-17c22/23o
GPALN_005799-T1	0	Y	n3-16c21/22o
GPALN_016181-T1	0	Y	n5-17c23/24o
GPALN_009730-T1	0	Y	n3-18c23/24o
GPALN_016298-T1	0	Y	n3-14c22/23o
GPALN_005901-T1	0	Y	n6-16c25/26o
GPALN_007196-T1	1	Y	n2-10c18/19o237-259i
GPALN_003793-T1	0	Y	n9-20c24/25o
GPALN_016098-T3	0	Y	n6-17c24/25o
GPALN_008791-T1	0	Y	n7-16c21/22o
GPALN_014133-T1	0	Y	n3-13c18/19o
GPALN_015531-T1	0	Y	n4-13c18/19o
GPALN_015075-T1	0	Y	n3-13c17/18o
GPALN_007921-T1	0	Y	n3-14c19/20o
GPALN_014713-T1	0	Y	n3-14c32/33o
GPALN_010787-T1	0	Y	n7-18c29/30o
GPALN_009822-T1	0	Y	n9-20c24/25o
GPALN_009571-T2	2	Y	n7-18c23/24o55-75i87-108o
GPALN_002919-T1	0	Y	n31-42c47/48o
GPALN_011807-T1	0	Y	n10-25c30/31o
GPALN_003379-T1	0	Y	n6-19c24/25o
GPALN_013301-T1	0	Y	n3-11c16/17o
GPALN_014962-T1	0	Y	n6-17c25/26o
GPALN_002133-T1	0	Y	n14-28c36/37o
GPALN_004930-T1	1	Y	n13-25c30/31o275-301i
GPALN_013805-T1	0	Y	n3-14c19/20o
GPALN_006766-T1	0	Y	n11-22c27/28o
GPALN_008458-T1	0	Y	n4-14c18/19o
GPALN_007368-T1	0	Y	n7-18c25/26o
GPALN_002201-T1	0	Y	n8-19c23/24o
GPALN_007345-T1	0	Y	n13-24c28/29o
GPALN_014523-T1	0	Y	n5-16c21/22o
GPALN_000870-T1	0	Y	n9-17c22/23o
GPALN_003092-T1	0	Y	n7-22c27/28o
GPALN_004376-T1	0	Y	n8-23c31/32o
GPALN_001572-T1	0	Y	n2-10c17/18o
GPALN_002746-T1	0	Y	n3-11c16/17o
GPALN_010416-T1	0	Y	n7-19c24/25o
GPALN_008869-T1	2	Y	n14-19c23/24o47-70i100-116o
GPALN_002168-T1	0	Y	n8-19c28/29o
GPALN_001144-T1	0	Y	n6-17c25/26o
GPALN_006778-T1	0	Y	n7-18c23/24o
GPALN_011435-T1	5	Y	n3-14c22/23o32-56i116-135o189-213i225-247o253-273i
GPALN_012285-T1	0	Y	n2-12c20/21o
GPALN_003860-T1	0	Y	n7-20c25/26o
GPALN_009925-T1	0	Y	n7-18c23/24o
GPALN_001760-T1	0	Y	n5-16c21/22o
GPALN_005770-T1	0	Y	n15-22c30/31o
GPALN_011611-T1	0	Y	n6-19c24/25o
GPALN_013555-T1	0	Y	n11-21c26/27o
GPALN_005100-T1	0	Y	n8-19c24/25o
GPALN_008272-T1	0	Y	n5-12c17/18o
GPALN_006727-T1	0	Y	n6-14c18/19o
GPALN_008762-T1	1	Y	n8-19c26/27o495-522i
GPALN_013210-T1	0	Y	n2-12c16/17o
GPALN_001735-T1	0	Y	n12-23c28/29o
GPALN_005611-T1	0	Y	n7-18c23/24o
GPALN_008501-T1	0	Y	n2-13c18/19o
GPALN_007222-T1	0	Y	n8-19c26/27o
GPALN_016055-T1	1	Y	n5-16c22/23o1203-1220i
GPALN_012445-T2	0	Y	n19-30c38/39o
GPALN_003301-T1	0	Y	n5-16c20/21o
GPALN_000555-T1	11	Y	n2-12c17/18o201-222i263-285o305-328i357-377o415-436i443-465o499-518i530-552o564-584i596-619o639-662i
GPALN_015025-T1	0	Y	n10-18c26/27o
GPALN_001110-T1	0	Y	n5-16c23/24o
GPALN_000342-T1	0	Y	n3-15c19/20o
GPALN_001298-T1	0	Y	n5-15c27/28o
GPALN_014040-T1	6	Y	n3-11c16/17o40-58i115-137o157-177i189-210o216-237i249-267o
GPALN_011804-T1	0	Y	n4-15c20/21o
GPALN_007699-T1	0	Y	n4-11c15/16o
GPALN_014077-T1	0	Y	n15-26c34/35o
GPALN_014890-T1	1	Y	n6-18c23/24o39-60i
GPALN_003143-T1	0	Y	n4-14c19/20o
GPALN_011884-T1	0	Y	n3-14c22/23o
GPALN_016193-T1	0	Y	n7-18c22/23o
GPALN_008111-T1	0	Y	n4-15c20/21o
GPALN_000350-T1	0	Y	n6-25c30/31o
GPALN_002489-T1	0	Y	n6-14c19/20o
GPALN_013295-T1	0	Y	n3-14c19/20o
GPALN_006932-T1	0	Y	n4-15c23/24o
GPALN_004968-T1	0	Y	n2-12c18/19o
GPALN_000435-T1	0	Y	n10-19c23/24o
GPALN_000153-T2	0	Y	n20-33c38/39o
GPALN_015272-T1	1	Y	n3-16c21/22o254-272i
GPALN_001823-T1	0	Y	n2-13c20/21o
GPALN_003358-T1	0	Y	n187-192c197/198o
GPALN_005909-T1	0	Y	n13-24c33/34o
GPALN_012765-T1	4	Y	n30-41c49/50o73-99i143-164o176-193i256-283o
GPALN_011650-T1	0	Y	n6-17c24/25o
GPALN_016128-T1	0	Y	n6-17c24/25o
GPALN_010411-T1	0	Y	n6-17c22/23o
GPALN_002527-T1	0	Y	n9-19c23/24o
GPALN_012285-T2	0	Y	n2-12c20/21o
GPALN_009380-T1	1	Y	n3-14c18/19o189-207i
GPALN_015499-T1	1	Y	n9-20c29/30o373-393i
GPALN_005263-T1	0	Y	n7-18c23/24o
GPALN_006801-T1	0	Y	n12-23c27/28o
GPALN_001352-T1	0	Y	n4-11c16/17o
GPALN_001546-T1	2	Y	n3-16c24/25o54-75i87-105o
GPALN_002495-T1	1	Y	n6-18c23/24o39-57i
GPALN_015307-T1	0	Y	n8-19c27/28o
GPALN_000924-T1	1	Y	n3-14c22/23o410-431i
GPALN_014568-T1	0	Y	n3-14c21/22o
GPALN_006343-T1	0	Y	n2-12c17/18o
GPALN_013480-T1	0	Y	n5-20c24/25o
GPALN_001024-T1	0	Y	n2-15c24/25o
GPALN_007606-T1	0	Y	n11-19c24/25o
GPALN_009583-T1	0	Y	n5-16c21/22o
GPALN_011722-T1	0	Y	n5-16c22/23o
GPALN_009015-T1	0	Y	n2-12c16/17o
GPALN_004559-T1	0	Y	n3-14c19/20o
GPALN_004816-T1	0	Y	n2-13c18/19o
GPALN_004738-T1	0	Y	n3-15c19/20o
GPALN_015175-T1	0	Y	n6-17c22/23o
GPALN_009901-T1	0	Y	n12-23c28/29o
GPALN_013375-T1	0	Y	n7-16c24/25o
GPALN_013190-T1	0	Y	n54-61c66/67o
GPALN_012140-T1	0	Y	n5-16c21/22o
GPALN_006613-T1	0	Y	n7-18c26/27o
GPALN_004209-T1	0	Y	n10-21c26/27o
GPALN_016115-T1	1	Y	n3-14c18/19o632-657i
GPALN_007089-T1	5	Y	n13-24c36/37o359-381i393-410o422-440i751-769o775-795i
GPALN_009281-T1	5	Y	n4-15c23/24o866-888i1430-1450o1462-1480i1492-1515o1521-1539i
GPALN_013712-T1	0	Y	n4-14c19/20o
GPALN_008337-T1	1	Y	n2-12c18/19o45-71i
GPALN_003941-T1	0	Y	n5-14c26/27o
GPALN_011385-T1	11	Y	n3-13c18/19o97-118i125-145o151-173i185-205o211-236i274-298o318-338i350-371o377-398i410-431o443-464i
GPALN_015640-T1	0	Y	n15-27c32/33o
GPALN_007850-T1	0	Y	n3-11c20/21o
GPALN_005732-T1	0	Y	n2-13c17/18o
GPALN_001745-T1	0	Y	n122-133c138/139o
GPALN_003534-T1	0	Y	n12-19c24/25o
GPALN_015521-T1	0	Y	n40-51c55/56o
GPALN_015181-T1	0	Y	n8-19c24/25o
GPALN_011627-T2	0	Y	n4-17c25/26o
GPALN_014171-T1	0	Y	n18-29c33/34o
GPALN_010137-T1	0	Y	n8-18c23/24o
GPALN_000526-T1	0	Y	n12-23c31/32o
GPALN_012366-T1	0	Y	n3-14c18/19o
GPALN_004951-T1	0	Y	n6-17c21/22o
GPALN_000449-T1	1	Y	n6-17c27/28o51-72i
GPALN_001496-T1	0	Y	n10-21c30/31o
GPALN_002177-T1	0	Y	n26-37c42/43o
GPALN_004480-T1	0	Y	n7-15c23/24o
GPALN_010583-T1	0	Y	n3-13c20/21o
GPALN_000704-T1	3	Y	n3-14c32/33o384-404i618-641o677-702i
GPALN_009029-T2	2	Y	n8-18c25/26o476-496i625-651o
GPALN_007488-T1	1	Y	n5-23c29/30o79-101i
GPALN_011304-T1	4	Y	n13-25c30/31o261-285i292-314o326-345i480-498o
GPALN_010425-T1	0	Y	n7-18c36/37o
GPALN_002348-T1	0	Y	n10-19c24/25o
GPALN_010629-T1	0	Y	n2-13c17/18o
GPALN_001806-T1	0	Y	n17-25c33/34o
GPALN_015326-T1	0	Y	n9-16c21/22o
GPALN_013252-T1	0	Y	n8-19c23/24o
GPALN_013639-T1	0	Y	n2-13c18/19o
GPALN_012116-T1	2	Y	n4-12c17/18o33-54i66-94o
GPALN_003951-T1	0	Y	n4-15c20/21o
GPALN_004800-T1	0	Y	n4-19c28/29o
GPALN_013438-T1	0	Y	n8-18c30/31o
GPALN_007386-T1	0	Y	n8-17c22/23o
GPALN_004857-T1	0	Y	n7-17c22/23o
GPALN_003697-T1	1	Y	n2-12c17/18o33-55i
GPALN_002552-T1	0	Y	n11-22c27/28o
GPALN_004739-T1	0	Y	n3-15c19/20o
GPALN_006223-T1	0	Y	n5-15c26/27o
GPALN_000179-T1	0	Y	n15-25c29/30o
GPALN_014341-T1	0	Y	n3-14c19/20o
GPALN_013552-T1	0	Y	n12-27c35/36o
GPALN_016040-T1	0	Y	n4-11c16/17o
GPALN_014398-T1	0	Y	n3-13c17/18o
GPALN_002750-T1	0	Y	n7-18c23/24o
GPALN_012452-T1	0	Y	n3-14c19/20o
GPALN_011364-T1	0	Y	n3-13c18/19o
GPALN_007468-T1	1	Y	n3-14c22/23o681-698i
GPALN_005181-T1	0	Y	n10-25c34/35o
GPALN_009682-T1	0	Y	n12-22c32/33o
GPALN_003447-T1	0	Y	n3-14c19/20o
GPALN_001131-T1	0	Y	n10-20c25/26o
GPALN_002299-T1	0	Y	n2-12c17/18o
GPALN_009912-T1	0	Y	n4-16c21/22o
GPALN_015711-T1	0	Y	n6-19c25/26o
GPALN_001770-T1	1	Y	n3-13c17/18o128-145i
GPALN_006381-T1	1	Y	n4-14c19/20o545-571i
GPALN_007739-T1	0	Y	n4-14c20/21o
GPALN_016168-T1	0	Y	n18-28c35/36o
GPALN_009660-T1	0	Y	n6-17c22/23o
GPALN_014897-T1	0	Y	n4-15c23/24o
GPALN_001186-T1	1	Y	n2-12c17/18o444-461i
GPALN_010970-T1	0	Y	n4-15c20/21o
GPALN_001738-T1	0	Y	n13-23c31/32o
GPALN_001356-T1	0	Y	n5-15c19/20o
GPALN_010297-T1	0	Y	n8-16c21/22o
GPALN_013575-T1	0	Y	n6-17c25/26o
GPALN_004604-T1	0	Y	n3-14c22/23o
GPALN_009358-T1	0	Y	n2-12c17/18o
GPALN_006780-T1	0	Y	n5-15c19/20o
GPALN_015100-T1	0	Y	n3-14c18/19o
GPALN_002267-T1	0	Y	n10-20c24/25o
GPALN_007221-T1	0	Y	n7-20c28/29o
GPALN_006884-T1	0	Y	n10-23c30/31o
GPALN_010346-T1	0	Y	n7-17c28/29o
GPALN_008972-T1	0	Y	n18-29c35/36o
GPALN_003908-T1	1	Y	n8-16c21/22o305-323i
GPALN_013819-T1	1	Y	n7-18c23/24o459-482i
GPALN_003757-T1	0	Y	n3-14c19/20o
GPALN_002948-T1	0	Y	n4-12c17/18o
GPALN_016384-T1	1	Y	n4-12c24/25o40-64i
GPALN_003368-T1	0	Y	n6-16c27/28o
GPALN_001845-T1	3	Y	n7-18c23/24o33-53i65-87o150-171i
GPALN_002052-T1	0	Y	n4-14c19/20o
GPALN_014002-T1	0	Y	n4-15c20/21o
GPALN_013504-T1	0	Y	n2-9c13/14o
GPALN_002848-T1	0	Y	n8-23c29/30o
GPALN_002604-T1	0	Y	n10-29c33/34o
GPALN_004440-T1	0	Y	n14-25c29/30o
GPALN_009815-T1	0	Y	n4-15c20/21o
GPALN_006102-T1	0	Y	n5-16c24/25o
GPALN_002482-T1	0	Y	n11-22c27/28o
GPALN_006886-T1	0	Y	n11-23c27/28o
GPALN_015771-T1	0	Y	n4-15c22/23o
GPALN_002345-T1	0	Y	n4-16c21/22o
GPALN_007028-T2	3	Y	n4-18c30/31o440-458i1040-1066o1086-1106i
GPALN_007059-T1	1	Y	n7-17c22/23o126-152i
GPALN_007134-T1	0	Y	n27-34c38/39o
GPALN_003061-T1	1	Y	n10-17c26/27o160-178i
GPALN_002760-T1	0	Y	n5-16c21/22o
GPALN_008113-T1	0	Y	n4-15c26/27o
GPALN_003091-T1	0	Y	n6-16c20/21o
GPALN_009109-T1	1	Y	n11-26c31/32o229-253i
GPALN_009550-T1	1	Y	n8-18c23/24o62-80i
GPALN_015061-T1	0	Y	n6-19c24/25o
GPALN_016209-T1	1	Y	n9-20c25/26o637-654i
GPALN_002493-T1	0	Y	n4-15c20/21o
GPALN_015372-T1	0	Y	n3-14c22/23o
GPALN_010168-T1	0	Y	n4-18c23/24o
GPALN_004554-T1	0	Y	n4-14c19/20o
GPALN_012407-T1	1	Y	n3-14c19/20o214-232i
GPALN_003527-T1	0	Y	n16-31c40/41o
GPALN_008930-T1	0	Y	n12-19c24/25o
GPALN_003850-T1	0	Y	n7-20c25/26o
GPALN_009314-T1	0	Y	n8-18c22/23o
GPALN_001013-T1	0	Y	n8-18c23/24o
GPALN_010899-T1	0	Y	n10-20c24/25o
GPALN_010816-T1	0	Y	n20-31c36/37o
GPALN_016293-T1	0	Y	n9-20c25/26o
GPALN_006568-T1	0	Y	n3-11c15/16o
GPALN_000858-T1	1	Y	n8-16c28/29o395-421i
GPALN_007974-T1	0	Y	n6-16c26/27o
GPALN_000092-T2	3	Y	n9-20c24/25o135-152i164-185o205-225i
GPALN_005018-T1	0	Y	n6-17c22/23o
GPALN_006034-T1	0	Y	n6-14c19/20o
GPALN_014924-T1	0	Y	n5-16c21/22o
GPALN_000149-T1	7	Y	n6-14c22/23o497-517i529-549o569-589i601-622o642-665i691-709o715-733i
GPALN_000629-T1	0	Y	n14-29c37/38o
GPALN_005123-T1	2	Y	n6-17c25/26o49-65i86-110o
GPALN_015013-T1	0	Y	n3-14c18/19o
GPALN_011179-T1	0	Y	n4-15c20/21o
GPALN_013464-T1	1	Y	n4-15c23/24o277-297i
GPALN_007627-T1	15	Y	n3-14c19/20o250-271i283-300o320-340i347-370o390-410i459-477o489-509i521-541o580-599i620-640o646-669i681-703o735-754i775-793o805-831i
GPALN_007711-T1	0	Y	n6-16c25/26o
GPALN_003925-T1	0	Y	n10-21c26/27o
GPALN_001899-T1	0	Y	n23-32c37/38o
GPALN_010131-T1	0	Y	n5-15c19/20o
GPALN_006216-T1	0	Y	n3-11c20/21o
GPALN_004645-T1	2	Y	n15-25c30/31o224-241i248-270o
GPALN_011649-T1	0	Y	n6-16c21/22o
GPALN_007545-T1	0	Y	n4-14c32/33o
GPALN_003405-T1	0	Y	n5-12c17/18o
GPALN_005823-T1	0	Y	n4-15c21/22o
GPALN_014368-T1	0	Y	n5-16c21/22o
GPALN_005574-T1	1	Y	n9-21c29/30o243-271i
GPALN_012301-T1	0	Y	n11-21c26/27o
GPALN_008135-T1	0	Y	n6-17c25/26o
GPALN_009613-T1	0	Y	n4-15c20/21o
GPALN_003795-T1	0	Y	n9-20c24/25o
GPALN_014168-T1	0	Y	n12-19c27/28o
GPALN_014842-T1	0	Y	n16-27c31/32o
GPALN_009084-T1	1	Y	n6-25c30/31o364-386i
GPALN_005090-T1	0	Y	n8-19c24/25o
GPALN_000327-T1	0	Y	n25-36c41/42o
GPALN_008744-T1	1	Y	n21-32c37/38o290-315i
GPALN_016284-T1	0	Y	n3-11c19/20o
GPALN_012283-T1	0	Y	n7-18c26/27o
GPALN_015319-T1	0	Y	n6-17c22/23o
GPALN_007537-T1	0	Y	n6-18c24/25o
GPALN_001230-T1	0	Y	n2-13c17/18o
GPALN_001134-T1	0	Y	n8-16c21/22o
GPALN_007580-T1	0	Y	n3-15c19/20o
GPALN_001223-T1	0	Y	n19-30c34/35o
GPALN_003422-T1	0	Y	n4-12c17/18o
GPALN_012343-T1	0	Y	n3-14c22/23o
GPALN_012025-T1	0	Y	n8-20c25/26o
GPALN_013562-T1	1	Y	n19-30c35/36o467-488i
GPALN_012637-T1	1	Y	n5-16c21/22o666-684i
GPALN_010730-T1	0	Y	n14-27c31/32o
GPALN_009091-T1	1	Y	n10-25c30/31o155-177i
GPALN_001953-T1	1	Y	n8-19c28/29o464-485i
GPALN_010702-T1	1	Y	n8-19c27/28o519-544i
GPALN_004342-T1	0	Y	n8-18c22/23o
GPALN_007704-T1	0	Y	n4-15c20/21o
GPALN_000850-T1	1	Y	n10-20c26/27o58-90i
GPALN_005242-T1	0	Y	n9-19c24/25o
GPALN_014569-T1	0	Y	n2-13c17/18o
GPALN_002130-T1	1	Y	n4-15c22/23o272-289i
GPALN_012841-T1	2	Y	n10-21c26/27o180-200i240-258o
GPALN_002246-T1	0	Y	n4-18c26/27o
GPALN_003792-T1	0	Y	n8-19c27/28o
GPALN_010789-T1	0	Y	n10-21c25/26o
GPALN_004556-T1	0	Y	n2-13c18/19o
GPALN_015318-T1	0	Y	n3-14c19/20o
GPALN_013083-T1	1	Y	n22-30c35/36o1072-1090i
GPALN_010561-T2	0	Y	n3-13c20/21o
GPALN_004743-T1	1	Y	n6-17c21/22o222-240i
GPALN_014952-T1	0	Y	n8-18c23/24o
GPALN_013349-T1	0	Y	n8-19c27/28o
GPALN_010083-T1	0	Y	n3-14c18/19o
GPALN_004736-T1	0	Y	n4-15c20/21o
GPALN_001202-T1	0	Y	n3-14c22/23o
GPALN_007073-T1	0	Y	n4-13c18/19o
GPALN_005088-T1	0	Y	n8-19c24/25o
GPALN_002316-T1	1	Y	n4-15c23/24o176-197i
GPALN_006925-T1	0	Y	n7-14c19/20o
GPALN_011985-T1	0	Y	n8-16c21/22o
GPALN_004782-T1	0	Y	n3-14c18/19o
GPALN_015632-T1	0	Y	n5-20c24/25o
GPALN_015072-T1	1	Y	n24-34c39/40o63-82i
GPALN_002143-T1	0	Y	n3-14c22/23o
GPALN_013702-T1	2	Y	n11-19c24/25o475-497i543-561o
GPALN_005571-T1	0	Y	n4-11c18/19o
GPALN_012096-T1	0	Y	n5-17c25/26o
GPALN_004493-T1	0	Y	n4-14c19/20o
GPALN_009850-T1	0	Y	n8-17c22/23o
GPALN_010112-T1	0	Y	n4-15c20/21o
GPALN_014814-T1	1	Y	n8-18c28/29o411-432i
GPALN_005598-T1	0	Y	n12-23c31/32o
GPALN_004182-T1	0	Y	n10-20c28/29o
GPALN_001018-T2	0	Y	n4-12c16/17o
GPALN_001206-T1	0	Y	n23-33c38/39o
GPALN_008761-T1	3	Y	n3-14c23/24o47-65i77-98o110-130i
GPALN_014552-T1	0	Y	n3-14c19/20o
GPALN_001145-T1	0	Y	n5-17c25/26o
GPALN_009628-T1	0	Y	n8-19c23/24o
GPALN_010798-T1	0	Y	n5-16c28/29o
GPALN_006602-T1	0	Y	n6-20c25/26o
GPALN_010093-T1	0	Y	n2-12c16/17o
GPALN_010519-T1	0	Y	n5-16c23/24o
GPALN_013248-T1	6	Y	n10-21c30/31o108-133i145-164o170-188i200-220o240-264i346-364o
GPALN_007195-T1	1	Y	n3-14c18/19o354-383i
GPALN_008216-T1	6	Y	n4-10c18/19o333-352i410-427o439-460i563-583o589-606i613-634o
GPALN_015926-T1	0	Y	n10-22c26/27o
GPALN_003381-T1	0	Y	n3-11c15/16o
GPALN_010170-T1	0	Y	n5-18c30/31o
GPALN_009585-T1	0	Y	n11-21c26/27o
GPALN_004000-T1	0	Y	n6-15c20/21o
GPALN_006300-T1	0	Y	n2-17c21/22o
GPALN_012887-T1	0	Y	n3-14c19/20o
GPALN_004292-T1	1	Y	n4-15c19/20o342-361i
GPALN_013415-T1	0	Y	n9-14c19/20o
GPALN_005647-T1	11	Y	n14-26c30/31o109-128i149-173o179-201i213-232o238-254i556-572o584-607i619-642o654-675i687-710o716-735i
GPALN_010535-T1	0	Y	n3-14c18/19o
GPALN_004346-T1	0	Y	n3-13c18/19o
GPALN_009910-T1	0	Y	n7-17c22/23o
GPALN_006425-T1	0	Y	n2-12c16/17o
GPALN_014727-T1	0	Y	n16-27c32/33o
GPALN_011411-T3	0	Y	n4-15c25/26o
GPALN_014903-T1	0	Y	n7-18c23/24o
GPALN_002692-T1	5	Y	n7-17c22/23o621-640i647-668o688-710i717-738o758-779i
GPALN_003956-T1	1	Y	n2-13c17/18o969-988i
GPALN_007556-T1	0	Y	n8-20c28/29o
GPALN_013793-T1	0	Y	n6-17c25/26o
GPALN_003794-T1	0	Y	n10-20c24/25o
GPALN_010636-T1	0	Y	n4-14c18/19o
GPALN_007317-T1	1	Y	n12-23c31/32o250-274i
GPALN_016347-T1	0	Y	n8-19c26/27o
GPALN_001846-T2	0	Y	n3-14c19/20o
GPALN_011209-T1	0	Y	n3-13c18/19o
GPALN_009325-T1	1	Y	n23-33c37/38o162-185i
GPALN_011867-T1	0	Y	n6-16c21/22o
GPALN_015364-T1	0	Y	n5-13c17/18o
GPALN_000274-T1	0	Y	n4-13c18/19o
GPALN_015636-T1	0	Y	n7-17c22/23o
GPALN_005175-T1	0	Y	n4-15c23/24o
GPALN_003819-T1	0	Y	n12-23c28/29o
GPALN_010307-T1	0	Y	n7-18c23/24o
GPALN_008335-T1	1	Y	n8-19c23/24o142-161i
GPALN_015602-T1	0	Y	n2-10c15/16o
GPALN_004166-T1	1	Y	n6-17c25/26o226-250i
GPALN_009848-T1	0	Y	n9-19c24/25o
GPALN_004380-T1	0	Y	n2-13c21/22o
GPALN_008681-T1	0	Y	n5-15c20/21o
GPALN_000362-T1	0	Y	n2-13c17/18o
GPALN_012829-T1	3	Y	n4-14c19/20o242-264i271-287o307-329i
GPALN_003970-T1	0	Y	n5-15c20/21o
GPALN_000196-T1	1	Y	n17-28c33/34o160-184i
GPALN_006417-T1	0	Y	n7-20c28/29o
GPALN_012667-T1	0	Y	n11-18c23/24o
GPALN_013482-T1	0	Y	n4-15c22/23o
GPALN_016017-T1	9	Y	n8-17c25/26o49-69i76-101o113-133i145-166o242-261i311-332o338-360i380-399o405-424i
GPALN_003222-T1	0	Y	n3-14c21/22o
GPALN_007142-T1	0	Y	n4-15c20/21o
GPALN_015236-T1	0	Y	n8-19c24/25o
GPALN_005087-T1	0	Y	n10-20c27/28o
GPALN_002221-T1	0	Y	n9-23c28/29o
GPALN_007605-T1	0	Y	n3-13c17/18o
GPALN_011854-T1	0	Y	n13-24c29/30o
GPALN_009596-T1	6	Y	n6-18c23/24o51-70i90-115o127-146i167-192o198-218i239-260o
GPALN_004254-T1	0	Y	n8-20c38/39o
GPALN_013269-T1	3	Y	n4-15c23/24o895-918i964-981o2105-2128i
GPALN_006963-T1	0	Y	n2-13c20/21o
GPALN_004464-T1	0	Y	n4-17c21/22o
GPALN_006988-T1	0	Y	n2-13c18/19o
GPALN_013720-T1	0	Y	n6-18c22/23o
GPALN_014957-T1	1	Y	n2-12c20/21o213-233i
GPALN_000556-T1	0	Y	n10-21c26/27o
GPALN_012979-T1	0	Y	n2-12c24/25o
GPALN_008771-T1	1	Y	n6-17c25/26o179-198i
GPALN_012544-T1	0	Y	n5-12c20/21o
GPALN_014044-T1	0	Y	n5-16c21/22o
GPALN_013637-T1	0	Y	n10-22c27/28o
GPALN_011027-T1	0	Y	n14-20c24/25o
GPALN_000308-T1	0	Y	n3-14c19/20o
GPALN_010388-T1	0	Y	n4-11c16/17o
GPALN_003714-T1	1	Y	n6-18c23/24o205-231i
GPALN_000692-T1	0	Y	n4-19c24/25o
GPALN_004986-T1	0	Y	n8-20c29/30o
GPALN_014843-T1	0	Y	n13-26c31/32o
GPALN_011814-T1	0	Y	n4-12c17/18o
GPALN_001862-T1	3	Y	n3-14c18/19o171-191i200-217o321-346i
GPALN_013109-T1	0	Y	n3-18c23/24o
GPALN_000163-T1	0	Y	n15-26c31/32o
GPALN_004626-T1	0	Y	n4-13c18/19o
GPALN_011571-T1	0	Y	n2-9c14/15o
GPALN_004840-T1	0	Y	n7-17c22/23o
GPALN_001153-T1	0	Y	n4-15c20/21o
GPALN_014929-T1	0	Y	n2-12c21/22o
GPALN_007185-T1	0	Y	n6-17c22/23o
GPALN_004051-T1	0	Y	n5-10c16/17o
GPALN_006384-T1	0	Y	n3-15c23/24o
GPALN_009513-T1	0	Y	n3-15c19/20o
GPALN_014061-T1	0	Y	n2-13c21/22o
GPALN_007574-T1	1	Y	n3-11c18/19o262-282i
GPALN_012981-T1	0	Y	n15-27c35/36o
GPALN_010933-T1	0	Y	n6-14c19/20o
GPALN_003408-T1	4	Y	n4-12c16/17o32-54i66-88o94-115i127-146o
GPALN_004840-T3	0	Y	n7-17c22/23o
GPALN_000607-T1	0	Y	n5-15c20/21o
GPALN_003468-T1	1	Y	n12-23c28/29o1182-1206i
GPALN_003100-T1	0	Y	n81-88c93/94o
GPALN_013230-T1	0	Y	n3-13c18/19o
GPALN_003684-T1	0	Y	n2-12c17/18o
GPALN_009039-T1	0	Y	n11-21c26/27o
GPALN_000176-T1	3	Y	n2-12c17/18o379-398i410-432o452-475i
GPALN_002678-T1	0	Y	n2-12c17/18o
GPALN_010621-T2	1	Y	n3-13c18/19o126-146i
GPALN_004308-T1	0	Y	n6-18c23/24o
GPALN_000566-T1	0	Y	n5-13c18/19o
GPALN_014794-T1	11	Y	n13-23c31/32o67-88i100-119o125-147i159-181o187-207i286-304o337-357i369-387o393-421i433-450o462-485i
GPALN_007123-T1	0	Y	n3-14c18/19o
GPALN_004373-T1	0	Y	n3-14c19/20o
GPALN_012718-T1	0	Y	n5-16c21/22o
GPALN_014147-T1	0	Y	n3-13c17/18o
GPALN_012207-T1	0	Y	n4-19c23/24o
GPALN_012186-T1	0	Y	n7-18c23/24o
GPALN_005198-T1	0	Y	n7-15c20/21o
GPALN_005731-T1	0	Y	n2-12c17/18o
GPALN_010402-T1	1	Y	n5-16c21/22o463-483i
GPALN_007512-T1	1	Y	n3-14c19/20o182-203i
GPALN_015601-T1	0	Y	n2-10c15/16o
GPALN_009568-T1	0	Y	n5-15c20/21o
GPALN_001106-T1	0	Y	n7-17c21/22o
GPALN_015295-T1	0	Y	n5-20c24/25o
GPALN_009584-T1	0	Y	n6-17c21/22o
GPALN_008112-T1	0	Y	n4-15c26/27o
GPALN_003906-T1	0	Y	n9-20c28/29o
GPALN_011620-T1	0	Y	n2-12c17/18o
GPALN_007801-T1	0	Y	n4-15c20/21o
GPALN_011054-T1	1	Y	n8-18c23/24o1518-1538i
GPALN_003396-T1	0	Y	n4-15c20/21o
GPALN_002377-T1	0	Y	n5-16c21/22o
GPALN_004845-T1	0	Y	n4-15c19/20o
GPALN_009335-T1	0	Y	n4-16c24/25o
GPALN_002519-T1	0	Y	n5-18c23/24o
GPALN_007388-T1	0	Y	n4-22c26/27o
GPALN_011865-T2	0	Y	n7-18c27/28o
GPALN_003246-T1	5	Y	n3-11c15/16o92-111i132-153o159-183i203-221o227-246i
GPALN_010270-T1	0	Y	n13-21c26/27o
GPALN_004646-T1	1	Y	n4-15c20/21o419-437i
GPALN_016299-T1	2	Y	n15-24c29/30o207-224i231-253o
GPALN_009622-T1	0	Y	n5-16c21/22o
GPALN_005829-T1	0	Y	n10-21c30/31o
GPALN_011891-T1	1	Y	n8-23c30/31o277-302i
GPALN_006748-T1	0	Y	n8-18c26/27o
GPALN_012189-T1	0	Y	n7-18c23/24o
GPALN_007374-T1	0	Y	n2-13c18/19o
GPALN_005081-T1	0	Y	n8-19c24/25o
GPALN_004064-T1	0	Y	n7-18c23/24o
GPALN_015250-T1	0	Y	n33-45c53/54o
GPALN_000092-T1	3	Y	n9-20c24/25o135-152i164-185o205-225i
GPALN_004602-T1	0	Y	n2-12c16/17o
GPALN_012526-T1	1	Y	n43-56c62/63o189-206i
GPALN_001063-T1	1	Y	n11-24c31/32o229-251i
GPALN_008815-T1	0	Y	n5-15c23/24o
GPALN_002264-T1	0	Y	n6-17c25/26o
GPALN_007862-T1	0	Y	n15-30c39/40o
GPALN_006318-T1	0	Y	n7-17c25/26o
GPALN_005097-T1	0	Y	n8-19c24/25o
GPALN_007332-T1	0	Y	n10-23c31/32o
GPALN_001748-T2	0	Y	n10-19c24/25o
GPALN_012760-T1	1	Y	n11-22c29/30o131-151i
GPALN_008798-T1	3	Y	n2-12c22/23o46-64i253-276o413-431i
GPALN_006765-T1	4	Y	n4-12c16/17o329-352i364-382o394-419i872-896o
GPALN_005753-T1	1	Y	n6-17c21/22o84-103i
GPALN_013816-T1	1	Y	n5-16c20/21o178-200i
GPALN_014746-T1	0	Y	n15-26c30/31o
GPALN_001897-T1	0	Y	n10-19c24/25o
GPALN_005105-T1	0	Y	n8-19c24/25o
GPALN_013280-T1	0	Y	n5-16c21/22o
GPALN_003852-T1	0	Y	n4-15c21/22o
GPALN_011220-T1	0	Y	n11-22c27/28o
GPALN_007616-T1	0	Y	n3-14c18/19o
GPALN_007198-T1	0	Y	n2-9c14/15o
GPALN_002905-T1	4	Y	n2-13c19/20o321-343i350-366o386-408i522-540o
GPALN_004265-T1	0	Y	n3-14c18/19o
GPALN_014501-T1	0	Y	n8-18c23/24o
GPALN_003420-T1	1	Y	n6-17c21/22o356-380i
GPALN_000957-T1	6	Y	n4-16c21/22o31-49i61-78o98-119i139-159o165-182i194-213o
GPALN_001101-T1	1	Y	n6-17c25/26o597-621i
GPALN_011070-T1	0	Y	n9-20c24/25o
GPALN_015823-T1	0	Y	n4-16c21/22o
GPALN_016213-T2	0	Y	n7-18c23/24o
GPALN_004078-T1	0	Y	n10-20c25/26o
GPALN_014537-T1	4	Y	n4-15c19/20o29-48i104-126o200-224i286-310o
GPALN_013054-T1	1	Y	n4-15c20/21o47-65i
GPALN_009609-T1	0	Y	n7-19c24/25o
GPALN_012358-T1	0	Y	n3-18c23/24o
GPALN_001763-T1	1	Y	n2-12c17/18o52-80i
GPALN_005697-T1	1	Y	n3-15c23/24o591-616i
GPALN_012592-T1	0	Y	n2-12c20/21o
GPALN_009038-T1	0	Y	n2-13c19/20o
GPALN_007547-T1	0	Y	n7-18c22/23o
GPALN_006907-T1	1	Y	n7-17c25/26o1005-1023i
GPALN_009978-T1	0	Y	n22-32c37/38o
GPALN_015302-T1	0	Y	n4-15c20/21o
GPALN_001207-T1	0	Y	n2-13c17/18o
GPALN_007936-T1	1	Y	n2-9c18/19o42-67i
GPALN_012921-T1	0	Y	n10-23c28/29o
GPALN_005918-T1	0	Y	n4-15c20/21o
GPALN_011881-T1	0	Y	n5-16c22/23o
GPALN_008102-T1	0	Y	n3-18c23/24o
GPALN_009532-T1	0	Y	n4-19c23/24o
GPALN_015476-T1	0	Y	n10-22c27/28o
GPALN_004398-T1	0	Y	n8-19c24/25o
GPALN_003186-T1	1	Y	n3-10c15/16o31-55i
GPALN_005538-T1	0	Y	n4-15c19/20o
GPALN_011319-T1	0	Y	n5-12c17/18o
GPALN_014005-T1	0	Y	n3-15c19/20o
GPALN_010299-T1	0	Y	n4-13c18/19o
GPALN_010581-T3	3	Y	n5-12c16/17o165-184i191-213o225-250i
GPALN_008089-T1	0	Y	n7-16c21/22o
GPALN_003182-T1	1	Y	n10-21c25/26o272-290i
GPALN_013282-T1	1	Y	n11-22c27/28o103-123i
GPALN_009872-T1	0	Y	n122-127c145/146o
GPALN_011023-T1	1	Y	n13-20c25/26o1105-1126i
GPALN_010568-T1	0	Y	n7-18c26/27o
GPALN_003415-T1	0	Y	n5-15c27/28o
GPALN_005867-T1	0	Y	n4-14c21/22o
GPALN_006594-T1	1	Y	n3-14c19/20o174-192i
GPALN_004707-T1	0	Y	n5-14c19/20o
GPALN_000680-T1	1	Y	n8-19c24/25o444-463i
GPALN_008816-T1	0	Y	n9-19c24/25o
GPALN_002450-T1	1	Y	n3-11c17/18o149-174i
GPALN_005527-T1	3	Y	n11-26c38/39o62-78i85-104o177-200i
GPALN_002762-T1	0	Y	n6-17c25/26o
GPALN_016178-T1	0	Y	n4-13c19/20o
GPALN_007173-T1	0	Y	n2-10c18/19o
GPALN_004482-T1	0	Y	n10-21c27/28o
GPALN_000352-T1	0	Y	n8-13c21/22o
GPALN_007115-T1	0	Y	n3-14c19/20o
GPALN_016265-T1	0	Y	n7-18c23/24o
GPALN_007905-T1	1	Y	n2-20c25/26o621-644i
GPALN_014271-T1	0	Y	n4-15c20/21o
GPALN_006629-T1	1	Y	n5-16c21/22o91-111i
GPALN_010475-T1	0	Y	n2-13c21/22o
GPALN_006328-T1	0	Y	n2-12c17/18o
GPALN_014182-T1	0	Y	n4-15c22/23o
GPALN_007048-T1	0	Y	n4-15c20/21o
GPALN_011797-T1	0	Y	n6-19c24/25o
GPALN_015365-T1	0	Y	n8-17c21/22o
GPALN_004992-T1	1	Y	n4-15c25/26o475-501i
GPALN_006843-T1	0	Y	n2-12c17/18o
GPALN_015234-T1	0	Y	n39-50c55/56o
GPALN_009791-T1	0	Y	n7-14c19/20o
GPALN_011030-T1	0	Y	n8-19c24/25o
GPALN_011143-T1	0	Y	n4-14c19/20o
GPALN_008658-T1	0	Y	n4-19c24/25o
GPALN_015299-T1	0	Y	n3-14c22/23o
GPALN_006061-T1	0	Y	n3-14c23/24o
GPALN_007411-T1	0	Y	n6-17c22/23o
GPALN_011439-T1	1	Y	n10-28c38/39o365-388i
GPALN_003837-T1	0	Y	n7-18c22/23o
GPALN_005061-T1	1	Y	n17-29c34/35o102-123i
GPALN_013937-T1	0	Y	n24-36c44/45o
GPALN_014075-T1	0	Y	n2-13c18/19o
GPALN_005160-T1	0	Y	n2-13c21/22o
GPALN_015736-T1	1	Y	n3-18c30/31o139-158i
GPALN_003218-T1	0	Y	n25-37c42/43o
GPALN_002084-T1	0	Y	n6-18c26/27o
GPALN_013490-T1	2	Y	n7-18c23/24o217-237i249-269o
GPALN_015188-T1	0	Y	n19-32c37/38o
GPALN_008502-T1	0	Y	n5-16c21/22o
GPALN_014503-T1	0	Y	n8-18c23/24o
GPALN_003120-T1	0	Y	n8-19c24/25o
GPALN_014498-T1	0	Y	n5-16c25/26o
GPALN_009441-T3	0	Y	n5-18c23/24o
GPALN_007592-T3	0	Y	n3-15c20/21o
GPALN_002666-T1	1	Y	n8-19c27/28o104-123i
GPALN_009366-T1	0	Y	n9-19c25/26o
GPALN_002554-T1	1	Y	n8-18c26/27o851-874i
GPALN_003596-T1	1	Y	n21-32c36/37o284-307i
GPALN_006471-T1	1	Y	n2-9c16/17o753-777i
GPALN_015174-T1	0	Y	n4-15c20/21o
GPALN_009146-T1	0	Y	n2-11c18/19o
GPALN_014530-T1	0	Y	n7-18c23/24o
GPALN_002398-T1	0	Y	n3-11c16/17o
GPALN_009441-T2	0	Y	n5-18c23/24o
GPALN_000576-T1	0	Y	n3-11c15/16o
GPALN_001094-T1	4	Y	n3-14c19/20o2272-2296i2361-2385o2405-2425i2437-2459o
GPALN_004790-T1	0	Y	n4-15c21/22o
GPALN_007376-T1	0	Y	n3-22c26/27o
GPALN_010536-T1	0	Y	n2-13c18/19o
GPALN_001643-T1	0	Y	n6-17c22/23o
GPALN_001087-T1	0	Y	n15-28c33/34o
GPALN_010735-T1	0	Y	n10-21c25/26o
GPALN_013538-T1	0	Y	n3-13c17/18o
GPALN_000432-T1	6	Y	n14-21c26/27o144-162i399-423o429-452i464-483o495-512i659-677o
GPALN_010726-T1	0	Y	n4-13c17/18o
GPALN_011726-T1	0	Y	n3-16c21/22o
GPALN_011426-T1	12	Y	n4-14c19/20o241-264i356-376o415-435i600-625o631-653i665-687o699-715i760-780o874-897i909-927o939-962i1100-1124o
GPALN_007883-T1	0	Y	n25-36c45/46o
GPALN_011766-T1	0	Y	n15-27c32/33o
GPALN_002079-T1	7	Y	n3-14c23/24o87-110i122-142o162-182i203-222o253-276i306-324o354-374i
GPALN_000112-T1	1	Y	n18-28c33/34o231-255i
GPALN_012577-T1	1	Y	n7-18c25/26o253-274i
GPALN_015301-T1	0	Y	n5-20c24/25o
GPALN_002122-T1	0	Y	n146-154c160/161o
GPALN_014191-T1	0	Y	n8-19c24/25o
GPALN_016188-T1	0	Y	n2-10c15/16o
GPALN_007384-T1	0	Y	n11-23c33/34o
GPALN_005069-T1	0	Y	n2-12c17/18o
GPALN_014772-T1	7	Y	n3-14c21/22o71-93i127-146o158-178i190-213o233-252i273-292o298-316i
GPALN_016268-T1	0	Y	n3-14c22/23o
GPALN_006949-T1	0	Y	n6-16c23/24o
GPALN_013108-T1	0	Y	n3-18c23/24o
GPALN_005244-T1	0	Y	n4-15c20/21o
GPALN_016153-T1	0	Y	n10-18c23/24o
GPALN_000115-T1	1	Y	n13-21c39/40o421-440i
GPALN_009573-T1	0	Y	n5-15c20/21o
GPALN_004368-T1	0	Y	n12-24c42/43o
GPALN_009505-T1	0	Y	n3-14c19/20o
GPALN_002671-T1	0	Y	n5-16c21/22o
GPALN_016349-T1	2	Y	n4-15c22/23o199-223i230-250o
GPALN_013828-T1	0	Y	n10-22c27/28o
GPALN_007361-T1	0	Y	n3-13c18/19o
GPALN_001989-T1	0	Y	n9-20c25/26o
GPALN_002464-T1	0	Y	n4-12c17/18o
GPALN_015832-T1	1	Y	n8-18c22/23o154-175i
GPALN_006728-T1	0	Y	n6-17c25/26o
GPALN_002713-T1	0	Y	n4-11c20/21o
GPALN_007592-T2	0	Y	n3-15c20/21o
GPALN_015826-T1	0	Y	n11-22c35/36o
GPALN_007082-T1	0	Y	n12-23c27/28o
GPALN_003118-T1	9	Y	n3-14c23/24o411-432i622-646o658-677i684-700o706-725i737-757o793-817i845-866o886-909i
GPALN_006455-T1	1	Y	n4-15c20/21o346-364i
GPALN_000889-T1	0	Y	n5-15c20/21o
GPALN_013404-T1	1	Y	n7-17c22/23o234-254i
GPALN_013143-T1	1	Y	n10-18c26/27o36-58i
GPALN_011720-T1	1	Y	n10-21c29/30o91-116i
GPALN_015211-T1	0	Y	n5-15c19/20o
GPALN_015600-T1	0	Y	n7-20c25/26o
GPALN_006603-T1	0	Y	n6-20c25/26o
GPALN_011737-T1	0	Y	n5-16c21/22o
GPALN_007072-T1	0	Y	n7-18c22/23o
GPALN_009251-T1	3	Y	n2-10c15/16o25-43i212-236o297-325i
GPALN_015559-T1	1	Y	n18-29c34/35o44-66i
GPALN_009649-T1	0	Y	n5-16c29/30o
GPALN_005916-T1	0	Y	n6-17c21/22o
GPALN_004194-T1	0	Y	n14-25c30/31o
GPALN_010617-T1	1	Y	n2-12c16/17o294-317i
GPALN_013630-T1	0	Y	n13-24c31/32o
GPALN_006752-T1	0	Y	n11-22c27/28o
GPALN_005401-T1	1	Y	n5-13c18/19o181-200i
GPALN_012842-T1	3	Y	n3-14c19/20o73-97i109-132o160-185i
GPALN_012021-T1	0	Y	n3-14c20/21o
GPALN_005070-T1	0	Y	n13-24c33/34o
GPALN_004561-T1	0	Y	n8-19c28/29o
GPALN_006596-T1	0	Y	n4-15c20/21o
GPALN_007038-T1	0	Y	n5-18c23/24o
GPALN_011926-T1	0	Y	n11-22c30/31o
GPALN_006655-T1	0	Y	n2-13c21/22o
GPALN_015367-T1	0	Y	n4-14c24/25o
GPALN_008394-T1	0	Y	n3-14c22/23o
GPALN_003649-T1	0	Y	n4-15c20/21o
GPALN_012854-T1	0	Y	n8-20c25/26o
GPALN_010968-T1	0	Y	n4-15c20/21o
GPALN_004902-T1	1	Y	n8-18c25/26o178-197i
GPALN_015244-T1	0	Y	n4-17c22/23o
GPALN_014884-T1	0	Y	n8-19c24/25o
GPALN_009632-T1	3	Y	n14-25c33/34o109-135i197-214o234-259i
GPALN_005930-T1	0	Y	n8-18c23/24o
GPALN_004655-T1	0	Y	n7-16c21/22o
GPALN_014935-T1	0	Y	n8-19c27/28o
GPALN_004253-T1	0	Y	n9-16c24/25o
GPALN_013827-T1	5	Y	n13-24c32/33o56-78i90-109o159-176i188-207o213-230i
GPALN_002492-T1	5	Y	n4-19c23/24o814-837i844-866o878-900i912-932o1037-1054i
GPALN_005177-T1	7	Y	n3-14c26/27o1198-1224i1236-1255o1267-1290i1311-1334o1360-1381i1393-1412o1432-1452i
GPALN_004418-T1	0	Y	n4-15c20/21o
GPALN_000537-T1	0	Y	n26-37c43/44o
GPALN_005734-T1	1	Y	n4-15c20/21o168-190i
GPALN_006532-T1	0	Y	n7-18c26/27o
GPALN_003480-T1	0	Y	n2-9c14/15o
GPALN_005969-T1	0	Y	n5-16c21/22o
GPALN_006598-T1	0	Y	n4-18c30/31o
GPALN_009869-T1	0	Y	n2-12c21/22o
GPALN_012964-T1	1	Y	n6-17c22/23o255-275i
GPALN_002224-T1	0	Y	n9-23c28/29o
GPALN_007663-T1	0	Y	n9-20c27/28o
GPALN_007396-T1	7	Y	n4-14c18/19o37-58i79-100o106-124i236-255o287-307i319-343o363-383i
GPALN_005291-T1	1	Y	n16-26c34/35o627-652i
GPALN_004879-T1	0	Y	n8-16c23/24o
GPALN_000361-T1	1	Y	n3-11c19/20o754-781i
GPALN_010903-T1	0	Y	n3-14c19/20o
GPALN_014885-T1	0	Y	n8-23c29/30o
GPALN_003949-T1	0	Y	n10-21c26/27o
GPALN_014958-T1	0	Y	n5-15c20/21o
GPALN_000517-T1	1	Y	n4-19c27/28o108-131i
GPALN_003942-T1	0	Y	n10-21c26/27o
GPALN_000842-T1	1	Y	n2-12c16/17o176-194i
GPALN_006901-T1	0	Y	n6-13c18/19o
GPALN_011065-T1	0	Y	n6-16c21/22o
GPALN_014378-T1	0	Y	n2-12c18/19o
GPALN_011496-T1	0	Y	n8-16c21/22o
GPALN_009829-T1	0	Y	n8-17c22/23o
GPALN_002500-T1	0	Y	n4-15c23/24o
GPALN_006029-T1	0	Y	n4-19c24/25o
GPALN_002942-T1	0	Y	n13-23c30/31o
GPALN_007680-T1	0	Y	n16-26c31/32o
GPALN_003376-T1	0	Y	n6-19c31/32o
GPALN_004939-T1	0	Y	n10-20c29/30o
GPALN_003834-T1	0	Y	n14-23c35/36o
GPALN_005147-T1	0	Y	n4-14c22/23o
GPALN_002515-T1	0	Y	n6-17c29/30o
GPALN_010444-T1	1	Y	n8-16c20/21o646-665i
GPALN_014764-T1	0	Y	n15-25c33/34o
GPALN_008867-T1	1	Y	n12-23c31/32o331-349i
GPALN_016283-T1	1	Y	n6-19c28/29o55-73i
GPALN_015436-T1	5	Y	n3-14c19/20o140-161i173-194o225-246i258-279o355-376i
GPALN_007061-T1	1	Y	n9-19c24/25o180-201i
GPALN_016278-T2	0	Y	n9-19c31/32o
GPALN_002086-T1	0	Y	n6-17c22/23o
GPALN_002782-T3	0	Y	n27-35c39/40o
GPALN_015450-T1	1	Y	n13-23c28/29o352-381i
GPALN_011139-T1	0	Y	n13-28c33/34o
GPALN_005612-T1	0	Y	n7-17c25/26o
GPALN_003176-T1	1	Y	n4-14c22/23o615-640i
GPALN_008535-T1	0	Y	n5-16c21/22o
GPALN_005237-T1	1	Y	n3-14c21/22o71-94i
GPALN_011615-T3	0	Y	n4-14c18/19o
GPALN_014172-T1	1	Y	n4-17c25/26o446-466i
GPALN_016378-T1	0	Y	n8-18c22/23o
GPALN_007536-T1	1	Y	n3-11c18/19o279-298i
GPALN_011510-T1	3	Y	n12-22c26/27o306-328i340-357o369-392i
GPALN_011917-T1	0	Y	n7-18c25/26o
GPALN_006945-T1	0	Y	n4-19c24/25o
GPALN_000469-T1	0	Y	n4-11c15/16o
GPALN_007498-T1	0	Y	n4-14c18/19o
GPALN_000265-T1	0	Y	n2-10c15/16o
GPALN_008057-T1	2	Y	n11-22c30/31o46-67i277-296o
GPALN_005738-T1	0	Y	n8-16c20/21o
GPALN_011159-T1	7	Y	n7-17c22/23o516-546i558-578o598-622i643-663o701-723i744-768o788-806i
GPALN_001139-T1	5	Y	n3-10c14/15o230-257i269-290o310-335i347-368o393-415i
GPALN_002466-T1	9	Y	n4-15c19/20o438-455i694-715o721-743i755-781o815-835i873-890o896-919i940-962o968-992i
GPALN_011808-T1	1	Y	n5-16c22/23o106-131i
GPALN_002669-T1	0	Y	n6-21c26/27o
GPALN_014052-T1	0	Y	n9-17c21/22o
GPALN_008681-T2	0	Y	n5-15c20/21o
GPALN_010569-T1	0	Y	n7-14c26/27o
GPALN_011399-T1	0	Y	n7-18c25/26o
GPALN_011777-T1	0	Y	n9-23c28/29o
GPALN_013896-T1	0	Y	n18-29c33/34o
GPALN_013687-T1	0	Y	n2-13c25/26o
GPALN_011556-T1	0	Y	n4-11c16/17o
GPALN_014053-T1	0	Y	n13-23c27/28o
GPALN_009893-T1	3	Y	n7-17c22/23o843-860i867-886o898-918i
GPALN_013829-T1	0	Y	n8-19c25/26o
GPALN_005911-T1	0	Y	n5-16c21/22o
GPALN_015372-T2	0	Y	n3-14c22/23o
GPALN_003918-T1	0	Y	n13-20c27/28o
GPALN_005924-T1	0	Y	n12-22c26/27o
GPALN_003889-T1	0	Y	n3-13c18/19o
GPALN_008839-T1	0	Y	n4-15c24/25o
GPALN_000894-T1	0	Y	n10-20c25/26o
GPALN_013617-T1	0	Y	n2-12c17/18o
GPALN_004571-T1	0	Y	n7-18c23/24o
GPALN_007962-T1	0	Y	n10-20c25/26o
GPALN_009593-T1	0	Y	n4-13c18/19o
GPALN_010391-T1	0	Y	n9-20c28/29o
GPALN_016270-T1	0	Y	n7-14c19/20o
GPALN_009502-T1	1	Y	n3-13c18/19o111-134i
GPALN_015415-T1	0	Y	n5-15c20/21o
GPALN_013114-T1	0	Y	n4-14c19/20o
GPALN_004410-T1	0	Y	n5-15c23/24o
GPALN_006201-T1	0	Y	n10-21c26/27o
GPALN_002361-T1	0	Y	n5-15c20/21o
GPALN_003011-T1	0	Y	n5-16c21/22o
GPALN_004750-T1	1	Y	n4-15c19/20o65-87i
GPALN_013204-T1	0	Y	n6-17c21/22o
GPALN_014126-T1	0	Y	n4-12c17/18o
GPALN_014835-T1	1	Y	n7-18c24/25o566-583i
GPALN_007259-T1	0	Y	n8-18c25/26o
GPALN_002433-T2	1	Y	n8-19c24/25o283-306i
GPALN_011411-T2	0	Y	n4-15c25/26o
GPALN_008941-T1	0	Y	n5-15c33/34o
GPALN_015297-T1	0	Y	n5-18c23/24o
GPALN_011660-T1	0	Y	n4-11c16/17o
GPALN_014309-T1	0	Y	n11-22c34/35o
GPALN_006023-T1	0	Y	n3-14c19/20o
GPALN_000113-T1	3	Y	n2-12c17/18o27-51i166-184o240-264i
GPALN_003178-T1	0	Y	n4-15c20/21o
GPALN_014959-T1	0	Y	n15-26c44/45o
GPALN_012300-T1	0	Y	n11-21c26/27o
GPALN_008557-T1	1	Y	n9-19c25/26o41-62i
GPALN_005226-T1	1	Y	n5-15c21/22o52-75i
GPALN_012465-T1	0	Y	n2-10c15/16o
GPALN_006772-T1	0	Y	n6-17c21/22o
GPALN_013900-T1	11	Y	n4-13c25/26o57-84i300-322o328-353i360-383o403-427i434-461o522-540i731-752o758-777i832-852o858-880i
GPALN_014060-T1	0	Y	n4-16c25/26o
GPALN_010580-T1	0	Y	n4-16c21/22o
GPALN_003094-T1	0	Y	n7-18c25/26o
GPALN_005172-T1	1	Y	n12-23c28/29o845-871i
GPALN_002354-T1	0	Y	n8-19c24/25o
GPALN_010321-T1	0	Y	n3-18c23/24o
GPALN_014392-T1	0	Y	n7-15c20/21o
GPALN_013210-T2	0	Y	n4-15c20/21o
GPALN_010587-T1	1	Y	n5-12c16/17o187-212i
GPALN_004587-T1	2	Y	n6-17c22/23o1582-1602i1713-1736o
GPALN_009944-T1	1	Y	n4-12c16/17o644-663i
GPALN_003890-T1	1	Y	n4-11c16/17o256-274i
GPALN_002192-T1	0	Y	n4-11c19/20o
GPALN_004345-T1	11	Y	n5-16c20/21o57-78i90-110o116-139i151-174o180-201i264-283o303-324i331-354o360-385i406-425o431-449i
GPALN_008764-T1	1	Y	n4-17c25/26o441-461i
GPALN_005758-T1	0	Y	n2-10c15/16o
GPALN_004072-T1	0	Y	n3-14c19/20o
GPALN_011371-T1	0	Y	n8-19c28/29o
GPALN_001590-T1	5	Y	n2-9c14/15o78-95i102-121o133-157i169-190o196-212i
GPALN_010945-T1	0	Y	n6-15c20/21o
GPALN_012731-T1	1	Y	n3-14c32/33o354-373i
GPALN_010185-T1	0	Y	n5-16c20/21o
GPALN_005889-T1	0	Y	n3-16c21/22o
GPALN_003217-T1	1	Y	n4-16c20/21o293-312i
GPALN_011519-T1	0	Y	n7-19c24/25o
GPALN_013448-T2	0	Y	n4-15c20/21o
GPALN_006720-T1	0	Y	n5-16c24/25o
GPALN_002425-T1	0	Y	n6-17c22/23o
GPALN_009498-T1	0	Y	n3-13c18/19o
GPALN_014253-T1	1	Y	n3-11c16/17o26-43i
GPALN_007533-T1	1	Y	n6-21c26/27o50-70i
GPALN_009650-T1	0	Y	n5-16c21/22o
GPALN_008739-T1	4	Y	n2-12c21/22o45-64i293-326o346-363i370-388o
GPALN_002036-T1	11	Y	n3-14c18/19o42-64i76-96o124-145i166-185o191-217i244-262o282-301i310-328o334-356i368-387o393-413i
GPALN_007845-T1	0	Y	n2-9c18/19o
GPALN_013170-T1	0	Y	n12-22c27/28o
GPALN_013941-T1	1	Y	n6-18c23/24o164-185i
GPALN_002714-T1	0	Y	n4-16c21/22o
GPALN_005038-T1	0	Y	n4-15c20/21o
GPALN_013296-T1	0	Y	n8-20c25/26o
GPALN_014868-T1	0	Y	n5-16c21/22o
GPALN_009582-T1	0	Y	n6-17c22/23o
GPALN_014802-T1	0	Y	n6-16c21/22o
GPALN_000495-T1	7	Y	n6-17c22/23o259-279i291-314o326-352i364-382o402-425i446-468o474-495i
GPALN_014966-T1	0	Y	n5-16c24/25o
GPALN_010419-T1	0	Y	n2-12c17/18o
GPALN_014753-T1	0	Y	n3-16c21/22o
GPALN_008161-T1	0	Y	n4-15c23/24o
GPALN_010420-T1	0	Y	n2-11c15/16o
GPALN_010355-T1	0	Y	n10-21c26/27o
GPALN_005915-T1	0	Y	n3-15c20/21o
GPALN_011375-T1	0	Y	n7-17c28/29o
GPALN_005928-T1	0	Y	n6-21c28/29o
GPALN_012644-T1	0	Y	n2-9c13/14o
GPALN_007990-T1	1	Y	n4-16c21/22o470-487i
GPALN_007270-T1	0	Y	n3-14c18/19o
GPALN_011500-T1	0	Y	n4-11c16/17o
GPALN_002914-T1	3	Y	n15-26c35/36o457-476i840-861o1244-1267i
GPALN_004842-T1	0	Y	n3-16c21/22o
GPALN_003060-T1	0	Y	n7-17c22/23o
GPALN_013561-T1	0	Y	n3-14c21/22o
GPALN_014091-T1	0	Y	n9-20c24/25o
GPALN_004355-T1	0	Y	n13-20c27/28o
GPALN_003873-T1	0	Y	n5-16c21/22o
GPALN_002386-T1	0	Y	n20-31c36/37o
GPALN_003816-T1	0	Y	n2-12c19/20o
GPALN_007806-T1	1	Y	n3-13c19/20o242-259i
GPALN_001796-T1	0	Y	n3-15c19/20o
GPALN_006346-T1	1	Y	n15-26c31/32o523-546i
GPALN_002174-T1	0	Y	n8-19c31/32o
GPALN_012232-T1	0	Y	n8-17c22/23o
GPALN_000778-T1	7	Y	n3-14c18/19o221-243i255-271o318-341i361-378o398-424i445-465o507-525i
GPALN_004647-T1	0	Y	n3-14c22/23o
GPALN_002422-T1	2	Y	n2-20c27/28o496-515i671-695o
GPALN_005254-T1	1	Y	n9-27c31/32o1862-1887i
GPALN_001120-T1	0	Y	n8-27c31/32o
GPALN_005451-T1	0	Y	n6-17c23/24o
GPALN_011300-T1	0	Y	n19-27c32/33o
GPALN_003855-T2	0	Y	n6-17c25/26o
GPALN_010650-T1	0	Y	n7-18c23/24o
GPALN_002881-T1	0	Y	n3-21c32/33o
GPALN_000330-T1	9	Y	n28-39c50/51o60-84i96-117o137-170i215-239o259-283i320-337o447-477i489-509o529-547i
GPALN_004624-T2	0	Y	n7-18c22/23o
GPALN_000702-T2	0	Y	n10-23c30/31o
GPALN_008510-T1	0	Y	n8-17c21/22o
GPALN_003910-T1	1	Y	n8-16c24/25o305-323i
GPALN_011172-T1	0	Y	n3-14c19/20o
GPALN_012291-T1	1	Y	n2-12c20/21o133-153i
GPALN_005725-T1	0	Y	n5-13c17/18o
GPALN_015446-T1	0	Y	n8-15c23/24o
GPALN_006919-T1	0	Y	n3-14c22/23o
GPALN_002955-T1	1	Y	n3-14c19/20o141-159i
GPALN_003441-T1	1	Y	n5-10c18/19o28-50i
GPALN_003471-T1	0	Y	n20-31c49/50o
GPALN_016192-T1	1	Y	n4-15c20/21o277-295i
GPALN_013428-T1	4	Y	n4-12c16/17o26-44i56-81o87-108i117-135o
GPALN_007705-T1	0	Y	n3-14c18/19o
GPALN_002313-T1	1	Y	n9-19c24/25o203-224i
GPALN_005133-T1	0	Y	n5-12c18/19o
GPALN_009902-T1	0	Y	n12-23c28/29o
GPALN_013133-T1	0	Y	n5-12c17/18o
GPALN_007237-T1	0	Y	n4-14c18/19o
GPALN_011359-T1	0	Y	n2-12c17/18o
GPALN_007214-T1	0	Y	n5-17c22/23o
GPALN_007489-T1	0	Y	n2-9c14/15o
GPALN_005089-T1	0	Y	n9-20c24/25o
GPALN_008081-T1	0	Y	n4-12c17/18o
GPALN_010795-T1	0	Y	n10-21c26/27o
GPALN_003883-T1	0	Y	n6-17c22/23o
GPALN_007293-T1	0	Y	n8-15c20/21o
GPALN_004569-T1	0	Y	n3-15c20/21o
GPALN_016387-T1	0	Y	n16-27c32/33o
GPALN_005789-T1	0	Y	n4-12c16/17o
GPALN_013786-T1	0	Y	n18-26c31/32o
GPALN_000621-T1	0	Y	n4-19c24/25o
GPALN_006100-T1	0	Y	n9-20c26/27o
GPALN_007532-T1	0	Y	n5-13c18/19o
GPALN_001214-T1	0	Y	n9-20c25/26o
GPALN_014939-T1	2	Y	n4-14c25/26o35-52i64-86o
GPALN_013913-T1	0	Y	n8-19c28/29o
GPALN_010991-T1	0	Y	n7-17c22/23o
GPALN_008479-T1	0	Y	n6-16c24/25o
GPALN_012418-T1	3	Y	n7-20c24/25o34-52i64-87o93-116i
GPALN_012906-T1	9	Y	n27-39c44/45o77-99i106-129o135-153i174-197o203-221i289-310o330-347i356-374o450-470i
GPALN_015037-T1	0	Y	n6-16c24/25o
GPALN_013323-T1	0	Y	n3-9c14/15o
GPALN_014766-T1	1	Y	n27-38c46/47o70-94i
GPALN_014640-T1	0	Y	n2-20c24/25o
GPALN_014472-T1	3	Y	n3-14c19/20o35-61i73-96o102-121i
GPALN_002352-T1	0	Y	n10-20c24/25o
GPALN_010613-T1	0	Y	n5-16c21/22o
GPALN_003846-T1	0	Y	n11-18c23/24o
GPALN_010851-T1	0	Y	n4-11c19/20o
GPALN_015605-T1	0	Y	n6-17c22/23o
GPALN_008105-T1	0	Y	n6-21c26/27o
GPALN_003825-T1	0	Y	n3-14c18/19o
GPALN_013284-T1	0	Y	n3-13c17/18o
GPALN_002885-T1	1	Y	n16-29c34/35o115-138i
GPALN_005042-T1	1	Y	n9-20c25/26o518-541i
GPALN_004098-T1	0	Y	n5-16c20/21o
GPALN_007198-T2	0	Y	n2-9c14/15o
GPALN_004056-T1	0	Y	n9-20c26/27o
GPALN_001252-T1	0	Y	n8-19c26/27o
GPALN_007693-T1	0	Y	n4-15c19/20o
GPALN_003117-T1	0	Y	n34-45c54/55o
GPALN_000440-T1	0	Y	n6-14c19/20o
GPALN_006654-T1	0	Y	n6-17c25/26o
GPALN_005537-T1	3	Y	n2-13c19/20o66-85i145-169o239-258i
GPALN_015157-T1	0	Y	n3-11c15/16o
GPALN_004011-T1	0	Y	n5-13c21/22o
GPALN_001179-T1	0	Y	n3-14c19/20o
GPALN_005383-T1	0	Y	n3-13c17/18o
GPALN_001190-T1	2	Y	n12-19c24/25o248-267i274-294o
GPALN_008597-T1	1	Y	n2-13c18/19o203-227i
GPALN_005109-T1	0	Y	n15-26c31/32o
GPALN_008084-T1	10	Y	n5-16c32/33o48-66i73-94o100-122i134-155o167-185i273-294o306-329i336-357o369-390i402-421o
GPALN_004687-T1	0	Y	n7-17c22/23o
GPALN_005766-T1	0	Y	n7-18c22/23o
GPALN_005180-T1	0	Y	n3-14c28/29o
GPALN_005801-T1	3	Y	n4-15c23/24o119-138i678-701o707-731i
GPALN_006820-T1	0	Y	n8-16c23/24o
GPALN_006726-T1	0	Y	n6-17c25/26o
GPALN_010415-T1	0	Y	n6-13c20/21o
GPALN_011464-T2	0	Y	n6-17c25/26o
GPALN_012215-T1	0	Y	n5-16c21/22o
GPALN_009717-T1	1	Y	n10-21c28/29o154-177i
GPALN_005820-T1	0	Y	n6-17c29/30o
GPALN_007644-T1	1	Y	n13-23c28/29o500-524i
GPALN_002989-T1	0	Y	n6-16c21/22o
GPALN_005893-T1	0	Y	n4-15c19/20o
GPALN_003449-T1	0	Y	n10-21c25/26o
GPALN_007538-T1	3	Y	n3-11c15/16o25-48i69-101o121-144i
GPALN_012287-T1	0	Y	n4-14c19/20o
GPALN_015161-T1	0	Y	n5-16c21/22o
GPALN_004195-T1	0	Y	n13-20c27/28o
GPALN_005326-T1	0	Y	n10-19c27/28o
GPALN_003309-T1	0	Y	n6-16c20/21o
GPALN_015314-T1	1	Y	n4-15c20/21o382-401i
GPALN_012260-T1	0	Y	n6-17c25/26o
GPALN_001534-T1	0	Y	n3-14c19/20o
GPALN_008886-T1	0	Y	n20-31c43/44o
GPALN_001530-T1	1	Y	n5-16c34/35o263-288i
GPALN_007831-T1	1	Y	n20-31c39/40o211-235i
GPALN_001923-T1	1	Y	n10-21c33/34o106-131i
GPALN_009041-T1	0	Y	n9-19c24/25o
GPALN_010487-T1	0	Y	n12-22c27/28o
GPALN_011812-T1	0	Y	n4-15c20/21o
GPALN_016152-T1	0	Y	n4-17c21/22o
GPALN_004178-T1	0	Y	n3-14c22/23o
GPALN_003179-T1	0	Y	n8-19c25/26o
GPALN_014539-T1	2	Y	n5-17c22/23o589-608i620-642o
GPALN_008107-T2	1	Y	n4-17c22/23o2863-2885i
GPALN_006313-T1	1	Y	n4-15c20/21o221-242i
GPALN_011609-T2	0	Y	n5-13c18/19o
GPALN_006659-T1	0	Y	n5-16c21/22o
GPALN_003769-T1	0	Y	n5-13c17/18o
GPALN_014393-T1	0	Y	n10-17c21/22o
GPALN_006631-T1	0	Y	n6-16c22/23o
GPALN_007683-T1	0	Y	n3-14c19/20o
GPALN_014865-T1	0	Y	n3-18c26/27o
GPALN_010136-T1	0	Y	n2-10c15/16o
GPALN_007866-T1	0	Y	n4-16c21/22o
GPALN_003743-T1	0	Y	n8-18c23/24o
GPALN_006599-T1	0	Y	n4-17c29/30o
GPALN_003015-T1	0	Y	n6-21c33/34o
GPALN_003826-T1	0	Y	n4-15c24/25o
GPALN_004515-T1	1	Y	n7-17c25/26o41-59i
GPALN_010732-T1	0	Y	n7-18c22/23o
GPALN_002418-T1	0	Y	n2-12c17/18o
GPALN_014037-T1	0	Y	n7-22c27/28o
GPALN_008908-T1	0	Y	n12-22c27/28o
GPALN_010524-T1	0	Y	n5-19c24/25o
GPALN_011865-T1	0	Y	n7-18c27/28o
GPALN_005750-T1	0	Y	n6-16c21/22o
GPALN_005868-T1	0	Y	n7-18c23/24o
GPALN_007077-T1	0	Y	n10-21c25/26o
GPALN_016091-T1	0	Y	n3-13c21/22o
GPALN_014001-T1	0	Y	n3-14c22/23o
GPALN_010316-T1	0	Y	n3-18c23/24o
GPALN_005416-T1	0	Y	n4-19c25/26o
GPALN_015704-T1	0	Y	n6-14c18/19o
GPALN_005236-T1	0	Y	n15-27c31/32o
GPALN_015739-T1	3	Y	n6-17c25/26o117-139i151-176o207-228i
GPALN_005439-T1	0	Y	n11-23c31/32o
GPALN_002701-T1	4	Y	n7-18c23/24o706-733i1167-1187o1290-1309i1634-1653o
GPALN_002550-T1	0	Y	n8-15c20/21o
GPALN_015224-T1	0	Y	n9-20c32/33o
GPALN_012445-T1	1	Y	n19-30c38/39o185-205i
GPALN_002427-T1	0	Y	n2-12c19/20o
GPALN_015780-T1	0	Y	n10-21c25/26o
GPALN_004306-T1	0	Y	n6-18c23/24o
GPALN_002117-T1	0	Y	n2-12c16/17o
GPALN_001501-T1	1	Y	n5-16c21/22o111-130i
GPALN_016167-T1	0	Y	n5-16c28/29o
GPALN_007560-T1	0	Y	n6-17c22/23o
GPALN_014377-T1	0	Y	n10-18c23/24o
GPALN_013306-T1	0	Y	n3-10c19/20o
GPALN_000527-T1	0	Y	n6-16c24/25o
GPALN_015182-T1	0	Y	n10-21c26/27o
GPALN_004117-T2	4	Y	n10-23c28/29o185-206i248-268o274-292i304-324o
GPALN_011546-T1	0	Y	n7-20c24/25o
GPALN_010799-T1	0	Y	n3-13c21/22o
GPALN_010810-T1	0	Y	n4-14c21/22o
GPALN_014080-T1	0	Y	n5-16c21/22o
GPALN_004555-T1	0	Y	n4-14c19/20o
GPALN_012363-T1	0	Y	n7-18c26/27o
GPALN_008866-T1	0	Y	n16-26c31/32o
GPALN_014477-T1	0	Y	n3-14c18/19o
GPALN_008701-T1	0	Y	n8-19c23/24o
GPALN_009253-T1	0	Y	n2-13c20/21o
GPALN_005185-T1	0	Y	n9-19c27/28o
GPALN_011583-T1	1	Y	n12-23c28/29o326-353i
GPALN_005299-T1	0	Y	n8-26c31/32o
GPALN_007074-T1	0	Y	n4-13c17/18o
GPALN_005049-T1	0	Y	n16-27c32/33o
GPALN_011018-T1	0	Y	n6-20c25/26o
GPALN_007129-T1	0	Y	n11-19c23/24o
GPALN_009008-T1	23	Y	n12-23c32/33o56-72i110-134o154-177i198-223o243-271i292-317o337-365i386-416o436-459i480-505o525-553i574-599o619-647i668-693o713-741i762-792o812-835i856-881o901-929i950-975o995-1022i1043-1059o1065-1084i
GPALN_004071-T1	0	Y	n8-17c21/22o
GPALN_012123-T1	0	Y	n3-14c22/23o
GPALN_012674-T1	0	Y	n10-21c26/27o
GPALN_011109-T1	0	Y	n3-13c21/22o
GPALN_006716-T1	0	Y	n3-17c22/23o
GPALN_004025-T1	0	Y	n5-13c18/19o
GPALN_005818-T1	0	Y	n3-11c16/17o
GPALN_006369-T1	0	Y	n5-16c20/21o
GPALN_007168-T1	0	Y	n3-11c15/16o
GPALN_004377-T1	0	Y	n3-14c19/20o
GPALN_012709-T1	1	Y	n3-14c19/20o29-45i
GPALN_007113-T1	0	Y	n4-15c23/24o
GPALN_008841-T1	0	Y	n3-14c20/21o
GPALN_007617-T1	0	Y	n6-16c21/22o
GPALN_005628-T1	0	Y	n10-18c22/23o
GPALN_006581-T1	0	Y	n4-15c20/21o
GPALN_009056-T1	0	Y	n25-40c45/46o
GPALN_015288-T1	1	Y	n4-14c19/20o247-268i
GPALN_002196-T1	0	Y	n8-19c31/32o
GPALN_007748-T1	0	Y	n8-19c24/25o
GPALN_011609-T1	0	Y	n5-13c18/19o
GPALN_009911-T1	0	Y	n7-17c22/23o
GPALN_002012-T1	1	Y	n3-13c18/19o34-54i
GPALN_004497-T1	0	Y	n10-22c26/27o
GPALN_009205-T1	0	Y	n6-19c24/25o
GPALN_015608-T1	0	Y	n3-11c17/18o
GPALN_016131-T1	0	Y	n3-14c19/20o
GPALN_008336-T1	1	Y	n8-23c29/30o163-187i
GPALN_000136-T1	0	Y	n9-19c26/27o
GPALN_002131-T1	0	Y	n6-17c22/23o
GPALN_003937-T1	0	Y	n3-13c18/19o
GPALN_012734-T1	0	Y	n3-18c23/24o
GPALN_001883-T1	0	Y	n14-24c29/30o
GPALN_014079-T1	0	Y	n9-17c21/22o
GPALN_000401-T1	6	Y	n3-10c19/20o275-296i317-336o356-376i413-437o479-500i507-531o
GPALN_004180-T1	0	Y	n6-16c21/22o
GPALN_011918-T2	0	Y	n10-22c30/31o
GPALN_007510-T1	1	Y	n3-14c19/20o181-202i
GPALN_002759-T1	0	Y	n9-20c28/29o
GPALN_015156-T1	0	Y	n7-18c25/26o
GPALN_005003-T1	0	Y	n5-16c20/21o
GPALN_006130-T1	0	Y	n3-10c18/19o
GPALN_002350-T1	1	Y	n3-12c17/18o269-289i
GPALN_011852-T1	0	Y	n31-44c49/50o
GPALN_008738-T1	0	Y	n9-20c27/28o
GPALN_007127-T1	0	Y	n6-17c21/22o
GPALN_015712-T1	0	Y	n6-19c24/25o
GPALN_003621-T1	5	Y	n6-16c21/22o400-421i572-593o609-628i635-661o690-720i
GPALN_006761-T1	0	Y	n5-13c18/19o
GPALN_006099-T1	0	Y	n26-37c43/44o
GPALN_001174-T1	0	Y	n5-17c22/23o
GPALN_016090-T1	0	Y	n8-19c24/25o
GPALN_001980-T1	2	Y	n3-11c19/20o151-172i312-333o
GPALN_004811-T1	0	Y	n5-13c18/19o
GPALN_001115-T1	0	Y	n6-17c28/29o
GPALN_013383-T1	0	Y	n4-15c20/21o
GPALN_006866-T1	0	Y	n9-20c25/26o
GPALN_014087-T1	0	Y	n13-24c31/32o
GPALN_002312-T1	1	Y	n9-19c24/25o206-227i
GPALN_003839-T1	0	Y	n7-17c22/23o
GPALN_004197-T1	0	Y	n4-15c22/23o
GPALN_007740-T1	1	Y	n11-21c29/30o169-192i
GPALN_005373-T1	0	Y	n5-15c19/20o
GPALN_000384-T2	1	Y	n19-29c34/35o348-371i
GPALN_009555-T2	0	Y	n6-17c25/26o
GPALN_003851-T1	1	Y	n3-14c19/20o234-254i
GPALN_009627-T1	0	Y	n6-17c22/23o
GPALN_010687-T1	7	Y	n7-15c24/25o79-101i138-159o165-184i196-220o240-259i280-297o303-321i
GPALN_014226-T1	0	Y	n7-18c23/24o
GPALN_015180-T1	0	Y	n4-14c22/23o
GPALN_014264-T1	1	Y	n8-18c25/26o237-255i
GPALN_004031-T1	0	Y	n5-16c21/22o
GPALN_013737-T1	1	Y	n6-18c23/24o862-882i
GPALN_011823-T1	0	Y	n4-15c20/21o
GPALN_013465-T2	0	Y	n5-14c19/20o
GPALN_005749-T1	1	Y	n5-17c22/23o175-197i
GPALN_004017-T1	0	Y	n6-18c23/24o
GPALN_015357-T1	0	Y	n6-17c25/26o
GPALN_015218-T1	0	Y	n7-17c22/23o
GPALN_013533-T1	0	Y	n4-22c27/28o
GPALN_010347-T1	0	Y	n3-11c16/17o
GPALN_008881-T1	0	Y	n8-18c23/24o
GPALN_005462-T1	0	Y	n4-15c20/21o
GPALN_010012-T1	0	Y	n18-28c35/36o
GPALN_012077-T1	0	Y	n3-11c19/20o
GPALN_005015-T1	1	Y	n18-29c34/35o142-159i
GPALN_005939-T1	0	Y	n5-15c20/21o
GPALN_007568-T1	1	Y	n10-21c29/30o629-655i
GPALN_010393-T1	1	Y	n3-11c19/20o246-276i
GPALN_005092-T1	1	Y	n9-18c22/23o509-538i
GPALN_008914-T1	0	Y	n8-18c26/27o
GPALN_002471-T1	0	Y	n8-18c23/24o
GPALN_007211-T1	5	Y	n2-13c23/24o97-116i137-158o164-188i208-226o232-250i
GPALN_002443-T1	1	Y	n10-20c24/25o34-55i
GPALN_001475-T1	1	Y	n5-16c21/22o113-131i
GPALN_005399-T1	0	Y	n12-20c28/29o
GPALN_009723-T1	1	Y	n10-21c28/29o154-177i
GPALN_007132-T1	0	Y	n11-21c25/26o
GPALN_010889-T1	1	Y	n11-22c27/28o62-88i
GPALN_014240-T1	0	Y	n6-13c18/19o
GPALN_012595-T1	1	Y	n9-20c28/29o926-950i
GPALN_007670-T1	0	Y	n6-17c22/23o
GPALN_003427-T1	0	Y	n3-13c21/22o
GPALN_005106-T1	0	Y	n2-12c16/17o
GPALN_014845-T1	0	Y	n3-15c20/21o
GPALN_014076-T1	0	Y	n10-18c25/26o
GPALN_005617-T1	0	Y	n12-23c28/29o
GPALN_004224-T1	0	Y	n4-12c16/17o
GPALN_015616-T1	0	Y	n3-12c17/18o
GPALN_003824-T1	0	Y	n6-17c21/22o
GPALN_006805-T1	2	Y	n10-21c25/26o258-278i290-309o
GPALN_009437-T1	0	Y	n8-18c26/27o
GPALN_009621-T1	0	Y	n3-14c18/19o
GPALN_011652-T1	0	Y	n3-14c22/23o
GPALN_004691-T1	0	Y	n9-17c21/22o
GPALN_014866-T1	0	Y	n5-16c21/22o
GPALN_010833-T1	0	Y	n4-19c28/29o
GPALN_015790-T1	0	Y	n13-20c28/29o
GPALN_007362-T1	0	Y	n15-25c30/31o
GPALN_015014-T1	0	Y	n41-51c58/59o
GPALN_003980-T1	7	Y	n6-17c22/23o63-84i96-115o121-141i153-179o185-205i217-234o240-263i
GPALN_010118-T1	0	Y	n6-17c24/25o
GPALN_008543-T1	0	Y	n8-18c25/26o
GPALN_003360-T1	1	Y	n6-19c24/25o48-68i
GPALN_008542-T1	4	Y	n7-17c22/23o235-257i269-288o300-322i519-539o
GPALN_007298-T1	0	Y	n2-10c16/17o
GPALN_009233-T1	0	Y	n3-14c24/25o
GPALN_011638-T1	0	Y	n8-19c23/24o
GPALN_013235-T1	0	Y	n8-20c32/33o
GPALN_012619-T1	3	Y	n7-17c26/27o239-258i270-288o300-323i
GPALN_001504-T1	0	Y	n16-27c32/33o
GPALN_003897-T1	0	Y	n5-17c22/23o
GPALN_008281-T1	0	Y	n5-13c20/21o
GPALN_000484-T1	0	Y	n4-15c33/34o
GPALN_004010-T1	0	Y	n3-14c19/20o
GPALN_006028-T1	1	Y	n3-13c19/20o250-272i
GPALN_009579-T1	0	Y	n7-18c24/25o
GPALN_012791-T1	0	Y	n8-19c24/25o
GPALN_008298-T1	0	Y	n19-27c34/35o
GPALN_012581-T1	1	Y	n18-25c43/44o155-177i
GPALN_003577-T1	1	Y	n26-37c42/43o1202-1227i
GPALN_007922-T1	1	Y	n12-23c28/29o359-380i
GPALN_000393-T1	0	Y	n4-14c19/20o
GPALN_014459-T1	0	Y	n15-26c34/35o
GPALN_012134-T1	1	Y	n5-16c21/22o238-258i
GPALN_005986-T1	0	Y	n14-25c30/31o
GPALN_015174-T2	0	Y	n4-15c20/21o
GPALN_011313-T1	0	Y	n25-36c41/42o
GPALN_012205-T1	0	Y	n4-14c18/19o
GPALN_013676-T1	0	Y	n3-10c15/16o
GPALN_000359-T1	0	Y	n10-21c25/26o
GPALN_015131-T1	0	Y	n2-13c21/22o
GPALN_013453-T1	0	Y	n19-26c44/45o
GPALN_002991-T1	0	Y	n6-17c22/23o
GPALN_009669-T1	0	Y	n9-19c27/28o
GPALN_006819-T1	0	Y	n8-18c22/23o
GPALN_010540-T1	0	Y	n6-16c23/24o
GPALN_014145-T1	0	Y	n4-16c24/25o
GPALN_002384-T1	3	Y	n9-20c25/26o131-151i163-183o203-224i
GPALN_015381-T1	0	Y	n6-17c22/23o
GPALN_014883-T1	0	Y	n5-16c24/25o
GPALN_007981-T1	0	Y	n8-18c25/26o
GPALN_015237-T1	0	Y	n8-19c24/25o
GPALN_010296-T1	0	Y	n5-16c21/22o
GPALN_014376-T1	0	Y	n3-10c15/16o
GPALN_001012-T1	0	Y	n7-18c24/25o
GPALN_011523-T1	0	Y	n4-15c23/24o
GPALN_012031-T1	3	Y	n10-21c26/27o61-84i105-124o166-184i
GPALN_005234-T1	1	Y	n4-14c19/20o1489-1510i
GPALN_002214-T1	0	Y	n4-13c18/19o
GPALN_007188-T1	2	Y	n6-15c20/21o261-280i300-322o
GPALN_002453-T1	0	Y	n11-22c27/28o
GPALN_001002-T1	0	Y	n8-19c23/24o
GPALN_010154-T1	0	Y	n3-13c21/22o
GPALN_006005-T1	0	Y	n10-20c25/26o
GPALN_010231-T1	0	Y	n4-15c20/21o
GPALN_001535-T1	1	Y	n2-12c20/21o134-156i
GPALN_000069-T1	0	Y	n19-29c34/35o
GPALN_006656-T1	0	Y	n3-14c19/20o
GPALN_014667-T1	7	Y	n6-17c22/23o50-75i87-106o112-134i146-166o178-197i209-230o236-254i
GPALN_004382-T1	0	Y	n3-14c19/20o
GPALN_000395-T1	0	Y	n2-13c18/19o
GPALN_001966-T1	0	Y	n15-25c31/32o
GPALN_014571-T1	1	Y	n4-23c28/29o124-147i
GPALN_007433-T1	0	Y	n2-13c20/21o
GPALN_011634-T1	0	Y	n3-13c21/22o
GPALN_009979-T1	1	Y	n11-30c38/39o321-345i
GPALN_009444-T1	0	Y	n2-12c17/18o
GPALN_007062-T1	0	Y	n4-15c19/20o
GPALN_005067-T1	0	Y	n8-19c24/25o
GPALN_011323-T1	0	Y	n3-13c18/19o
GPALN_014971-T1	0	Y	n7-18c26/27o
GPALN_002165-T1	0	Y	n4-15c23/24o
GPALN_000429-T1	0	Y	n3-14c19/20o
GPALN_002294-T1	0	Y	n5-20c24/25o
GPALN_014669-T1	7	Y	n6-17c22/23o74-94i106-125o131-154i166-186o198-217i238-259o265-283i
GPALN_013384-T1	0	Y	n4-15c20/21o
GPALN_000017-T1	0	Y	n2-12c17/18o
GPALN_015708-T1	1	Y	n7-15c20/21o147-168i
GPALN_014773-T1	0	Y	n8-19c24/25o
GPALN_007277-T1	0	Y	n3-14c19/20o
GPALN_004091-T1	0	Y	n3-21c26/27o
GPALN_014051-T1	1	Y	n4-17c21/22o37-66i
GPALN_015815-T1	0	Y	n4-12c17/18o
GPALN_016157-T1	0	Y	n10-18c23/24o
GPALN_002564-T1	0	Y	n3-10c15/16o
GPALN_006839-T1	0	Y	n4-15c20/21o
GPALN_013509-T1	0	Y	n17-24c29/30o
GPALN_008878-T1	0	Y	n6-17c25/26o
GPALN_011636-T1	0	Y	n3-14c18/19o
GPALN_007514-T1	0	Y	n6-15c19/20o
GPALN_002148-T1	0	Y	n14-28c36/37o
GPALN_011033-T1	0	Y	n18-37c42/43o
GPALN_009589-T1	0	Y	n5-16c24/25o
GPALN_002176-T1	0	Y	n5-16c21/22o
GPALN_012183-T1	0	Y	n7-16c21/22o
GPALN_001647-T1	0	Y	n3-9c14/15o
GPALN_009333-T4	3	Y	n6-15c20/21o56-77i183-204o216-242i
GPALN_007217-T1	0	Y	n4-15c23/24o
GPALN_000718-T1	4	Y	n2-13c19/20o301-323i335-352o364-387i504-521o
GPALN_003946-T1	0	Y	n10-21c26/27o
GPALN_008323-T1	1	Y	n15-26c30/31o203-227i
GPALN_007732-T1	1	Y	n3-14c19/20o90-110i
GPALN_000707-T1	0	Y	n7-18c26/27o
GPALN_003952-T1	0	Y	n4-15c20/21o
GPALN_008801-T1	7	Y	n3-10c22/23o46-78i99-119o139-159i179-199o267-292i324-341o347-370i
GPALN_011309-T1	0	Y	n9-20c27/28o
GPALN_011005-T1	1	Y	n16-27c32/33o458-477i
GPALN_004858-T1	0	Y	n7-18c23/24o
GPALN_016368-T1	0	Y	n7-15c19/20o
GPALN_015425-T1	0	Y	n3-14c18/19o
GPALN_010414-T1	0	Y	n7-19c24/25o
GPALN_015542-T1	1	Y	n8-17c22/23o144-165i
GPALN_014804-T1	2	Y	n2-13c18/19o220-243i654-675o
GPALN_008695-T1	0	Y	n13-24c29/30o
GPALN_006770-T1	0	Y	n13-21c29/30o
GPALN_008553-T1	2	Y	n5-16c22/23o464-487i494-516o
GPALN_014047-T1	7	Y	n2-9c17/18o27-50i71-91o111-139i159-179o216-240i261-290o296-321i
GPALN_012838-T1	0	Y	n6-16c21/22o
GPALN_003957-T1	1	Y	n4-15c23/24o267-285i
GPALN_015407-T1	0	Y	n4-14c19/20o
GPALN_014068-T1	0	Y	n4-15c19/20o
GPALN_008108-T1	0	Y	n3-16c21/22o
GPALN_014814-T2	1	Y	n8-18c28/29o411-432i
GPALN_005315-T1	0	Y	n13-24c28/29o
GPALN_002080-T1	0	Y	n4-14c19/20o
GPALN_007838-T1	0	Y	n2-12c20/21o
GPALN_006027-T1	1	Y	n3-13c19/20o242-262i
GPALN_001121-T1	1	Y	n2-13c20/21o503-526i
GPALN_013528-T1	0	Y	n5-15c19/20o
GPALN_010074-T1	0	Y	n8-19c23/24o
GPALN_008050-T1	1	Y	n11-21c25/26o514-540i
GPALN_000035-T1	0	Y	n2-12c17/18o
GPALN_014357-T1	0	Y	n10-21c25/26o
GPALN_002864-T1	0	Y	n5-16c21/22o
GPALN_003071-T1	0	Y	n2-13c25/26o
GPALN_007643-T1	0	Y	n5-16c21/22o
GPALN_010836-T1	2	Y	n8-18c26/27o185-206i213-231o
GPALN_006445-T1	0	Y	n4-16c24/25o
GPALN_014304-T1	0	Y	n8-23c28/29o
GPALN_007796-T1	0	Y	n4-19c24/25o
GPALN_004316-T1	2	Y	n13-24c29/30o228-251i263-284o
GPALN_002908-T1	0	Y	n14-25c30/31o
GPALN_006063-T1	0	Y	n5-15c23/24o
GPALN_004699-T1	0	Y	n7-18c26/27o
GPALN_003510-T1	0	Y	n17-28c33/34o
GPALN_002222-T1	0	Y	n3-21c28/29o
GPALN_012489-T1	0	Y	n17-28c38/39o
GPALN_006812-T1	0	Y	n8-16c23/24o
GPALN_000106-T1	1	Y	n2-12c17/18o140-163i
GPALN_002494-T1	0	Y	n10-21c33/34o
GPALN_009675-T1	0	Y	n4-14c19/20o
GPALN_001796-T2	1	Y	n3-15c19/20o387-407i
GPALN_014146-T1	0	Y	n4-16c24/25o
GPALN_003730-T1	0	Y	n3-11c19/20o
GPALN_006775-T1	0	Y	n9-20c24/25o
GPALN_015361-T1	0	Y	n2-13c17/18o
GPALN_008401-T1	5	Y	n27-38c43/44o652-670i682-704o724-747i782-809o829-849i
GPALN_000360-T1	0	Y	n6-16c21/22o
GPALN_009879-T1	0	Y	n5-14c19/20o
GPALN_013115-T1	1	Y	n11-21c26/27o173-196i
GPALN_015655-T1	0	Y	n12-22c27/28o
GPALN_004371-T1	0	Y	n3-13c18/19o
GPALN_005174-T1	0	Y	n6-17c22/23o
GPALN_006482-T1	0	Y	n10-21c29/30o
GPALN_008426-T1	1	Y	n2-7c12/13o46-68i
GPALN_009458-T1	0	Y	n4-15c20/21o
GPALN_013107-T1	0	Y	n3-18c23/24o
GPALN_013158-T1	0	Y	n6-17c22/23o
GPALN_007619-T1	3	Y	n6-16c21/22o37-57i64-86o106-126i
GPALN_015193-T1	0	Y	n4-15c20/21o
GPALN_013465-T1	0	Y	n5-14c19/20o
GPALN_013580-T1	0	Y	n8-19c24/25o
GPALN_010408-T1	0	Y	n8-17c22/23o
GPALN_014372-T1	0	Y	n5-16c21/22o
GPALN_015654-T1	0	Y	n10-21c26/27o
GPALN_010298-T1	0	Y	n5-16c21/22o
GPALN_015123-T1	0	Y	n7-18c26/27o
GPALN_009510-T1	0	Y	n12-23c32/33o
GPALN_011069-T2	0	Y	n16-25c33/34o
GPALN_001780-T1	0	Y	n2-12c20/21o
GPALN_002380-T1	0	Y	n2-12c20/21o
GPALN_011693-T1	1	Y	n4-12c18/19o62-80i
GPALN_009576-T1	0	Y	n10-21c30/31o
GPALN_000790-T1	1	Y	n2-11c23/24o451-478i
GPALN_004013-T1	1	Y	n4-15c27/28o369-389i
GPALN_003959-T1	0	Y	n11-18c25/26o
GPALN_010166-T1	0	Y	n4-14c18/19o
GPALN_004878-T1	0	Y	n7-18c26/27o
GPALN_015222-T1	0	Y	n8-19c28/29o
GPALN_003144-T1	1	Y	n8-16c24/25o366-388i
GPALN_006853-T1	0	Y	n4-15c20/21o
GPALN_002547-T1	0	Y	n10-17c23/24o
GPALN_014329-T1	0	Y	n6-14c19/20o
GPALN_005491-T1	7	Y	n14-24c32/33o135-160i172-196o211-231i252-277o307-325i929-952o972-990i
GPALN_006319-T1	0	Y	n2-12c20/21o
GPALN_011548-T1	0	Y	n16-26c31/32o
GPALN_008260-T1	7	Y	n5-12c17/18o58-83i95-118o138-161i173-192o242-262i294-310o330-352i
GPALN_004668-T1	0	Y	n3-13c20/21o
GPALN_009663-T1	2	Y	n7-18c23/24o217-237i249-269o
GPALN_014038-T1	0	Y	n7-18c23/24o
GPALN_013660-T1	1	Y	n10-21c26/27o249-267i
GPALN_014352-T1	1	Y	n3-13c19/20o247-269i
GPALN_015012-T1	0	Y	n7-18c26/27o
GPALN_006053-T1	0	Y	n6-14c23/24o
GPALN_012318-T1	0	Y	n14-21c26/27o
GPALN_003073-T1	0	Y	n7-20c25/26o
GPALN_010520-T1	0	Y	n5-18c26/27o
GPALN_013821-T1	1	Y	n8-19c24/25o492-513i
GPALN_006575-T1	0	Y	n8-19c24/25o
GPALN_011394-T1	5	Y	n8-19c27/28o123-146i175-197o468-493i514-541o879-904i
GPALN_004350-T1	0	Y	n7-18c23/24o
GPALN_007422-T1	0	Y	n7-18c26/27o
GPALN_007684-T1	0	Y	n10-17c21/22o
GPALN_005107-T1	0	Y	n8-19c24/25o
GPALN_010648-T1	0	Y	n8-19c26/27o
GPALN_008366-T1	0	Y	n6-16c23/24o
GPALN_002880-T1	0	Y	n4-15c23/24o
GPALN_001628-T1	0	Y	n11-21c25/26o
GPALN_008855-T1	1	Y	n16-27c33/34o43-69i
GPALN_010038-T1	2	Y	n3-14c26/27o50-70i323-345o
GPALN_009225-T1	0	Y	n12-22c30/31o
GPALN_002223-T1	0	Y	n3-21c29/30o
GPALN_015073-T1	0	Y	n4-16c24/25o
GPALN_012435-T1	1	Y	n6-17c22/23o99-118i
GPALN_003580-T1	4	Y	n10-23c28/29o329-348i360-378o390-416i428-445o
GPALN_009748-T1	0	Y	n10-21c26/27o
GPALN_008846-T1	0	Y	n5-12c17/18o
GPALN_014351-T1	1	Y	n3-13c19/20o247-269i
GPALN_015938-T1	0	Y	n13-21c29/30o
GPALN_002366-T1	0	Y	n6-17c29/30o
GPALN_003711-T1	0	Y	n38-49c61/62o
GPALN_005711-T1	1	Y	n23-33c38/39o104-129i
GPALN_002600-T1	0	Y	n5-16c21/22o
GPALN_000518-T1	0	Y	n10-21c26/27o
GPALN_005872-T1	0	Y	n4-15c20/21o
GPALN_010793-T1	0	Y	n10-21c25/26o
GPALN_001881-T1	0	Y	n2-13c17/18o
GPALN_001644-T1	1	Y	n6-13c17/18o813-832i
GPALN_009074-T1	3	Y	n8-19c23/24o607-625i678-695o861-883i
GPALN_007215-T1	0	Y	n7-18c27/28o
GPALN_001515-T1	0	Y	n15-25c30/31o
GPALN_015358-T1	0	Y	n8-19c25/26o
GPALN_004794-T1	0	Y	n6-17c25/26o
GPALN_002782-T1	1	Y	n27-35c39/40o436-460i
GPALN_010458-T1	2	Y	n4-16c26/27o313-333i744-765o
GPALN_013272-T1	0	Y	n3-14c19/20o
GPALN_007622-T1	0	Y	n3-13c21/22o
GPALN_011188-T1	0	Y	n7-18c23/24o
GPALN_013615-T1	1	Y	n2-12c17/18o65-84i
GPALN_012989-T1	0	Y	n6-17c25/26o
GPALN_007286-T2	0	Y	n5-16c20/21o
GPALN_003760-T1	7	Y	n3-14c20/21o192-211i223-245o276-298i318-343o370-394i477-497o509-534i
GPALN_009374-T1	1	Y	n5-16c24/25o705-730i
GPALN_005189-T1	0	Y	n6-13c17/18o
GPALN_014938-T1	0	Y	n10-21c26/27o
GPALN_009892-T1	0	Y	n4-22c27/28o
GPALN_004167-T1	1	Y	n3-13c19/20o127-150i
GPALN_006663-T1	0	Y	n7-19c24/25o
GPALN_013254-T1	0	Y	n4-14c19/20o
GPALN_006030-T1	0	Y	n5-16c24/25o
GPALN_002288-T1	0	Y	n5-20c24/25o
GPALN_000514-T1	0	Y	n6-14c19/20o
GPALN_005536-T1	7	Y	n4-15c20/21o30-48i60-78o84-104i111-130o142-165i177-195o270-294i
GPALN_001355-T1	0	Y	n14-25c32/33o
GPALN_015687-T1	0	Y	n8-20c25/26o
GPALN_003308-T1	0	Y	n3-15c19/20o
GPALN_000558-T1	0	Y	n5-14c19/20o
GPALN_010073-T1	0	Y	n6-17c25/26o
GPALN_011918-T1	0	Y	n10-22c30/31o
GPALN_006687-T1	0	Y	n8-15c19/20o
GPALN_008887-T1	1	Y	n23-38c43/44o1861-1883i
GPALN_013385-T1	0	Y	n5-18c23/24o
GPALN_008228-T1	1	Y	n7-18c26/27o179-201i
GPALN_013144-T1	0	Y	n9-20c28/29o
GPALN_008861-T1	4	Y	n7-12c17/18o33-61i82-101o121-145i157-178o
GPALN_010884-T1	5	Y	n3-10c15/16o68-89i101-118o124-139i146-165o177-197i
GPALN_010737-T1	0	Y	n7-18c23/24o
GPALN_002960-T1	0	Y	n6-14c18/19o
GPALN_013884-T1	5	Y	n2-12c17/18o27-46i58-80o196-218i225-245o265-291i
GPALN_003632-T1	1	Y	n23-34c38/39o864-885i
GPALN_004689-T1	0	Y	n3-14c19/20o
GPALN_010943-T1	0	Y	n4-12c17/18o
GPALN_011377-T1	0	Y	n5-16c23/24o
GPALN_004447-T1	0	Y	n5-13c22/23o
GPALN_007404-T1	0	Y	n11-18c23/24o
GPALN_014649-T1	1	Y	n3-18c23/24o733-752i
GPALN_011438-T1	0	Y	n11-22c27/28o
GPALN_012998-T1	0	Y	n6-16c23/24o
GPALN_011343-T1	1	Y	n23-34c39/40o729-750i
GPALN_001338-T1	1	Y	n4-19c24/25o291-308i
GPALN_004810-T1	0	Y	n2-13c18/19o
GPALN_014261-T1	0	Y	n4-12c16/17o
GPALN_000740-T1	0	Y	n19-30c37/38o
GPALN_011955-T1	0	Y	n13-20c27/28o
GPALN_007641-T1	9	Y	n10-21c25/26o326-343i363-385o424-446i458-486o492-516i537-570o590-612i633-656o662-685i
GPALN_003223-T1	0	Y	n3-14c21/22o
GPALN_002330-T1	0	Y	n3-13c17/18o
GPALN_004582-T1	0	Y	n15-27c32/33o
GPALN_012457-T1	0	Y	n2-13c18/19o
GPALN_002257-T1	0	Y	n5-16c21/22o
GPALN_002884-T1	0	Y	n8-20c28/29o
GPALN_009273-T1	0	Y	n21-32c37/38o
GPALN_015622-T1	0	Y	n2-13c20/21o
GPALN_008110-T1	0	Y	n4-15c20/21o
GPALN_008860-T1	0	Y	n9-20c30/31o
GPALN_003228-T1	0	Y	n23-33c38/39o
GPALN_011545-T1	0	Y	n3-14c19/20o
GPALN_008106-T1	0	Y	n2-12c17/18o
GPALN_007557-T1	0	Y	n3-14c19/20o
GPALN_010255-T1	0	Y	n2-11c16/17o
GPALN_003954-T1	0	Y	n4-15c20/21o
GPALN_013496-T1	0	Y	n11-21c29/30o
GPALN_009292-T1	0	Y	n8-19c24/25o
GPALN_012254-T1	3	Y	n8-19c27/28o37-56i68-87o99-119i
GPALN_015484-T1	9	Y	n18-28c40/41o1161-1181i1334-1358o1364-1397i1418-1439o1451-1467i1487-1505o1592-1612i1753-1774o1794-1815i
GPALN_009903-T1	0	Y	n2-12c20/21o
GPALN_010853-T1	0	Y	n11-19c27/28o
GPALN_010598-T1	0	Y	n3-14c32/33o
GPALN_013005-T1	1	Y	n6-17c25/26o102-121i
GPALN_005424-T1	0	Y	n17-28c32/33o
GPALN_007023-T1	0	Y	n3-11c29/30o
GPALN_006401-T1	0	Y	n11-21c33/34o
GPALN_003416-T2	0	Y	n5-13c18/19o
GPALN_006589-T1	0	Y	n6-17c23/24o
GPALN_004124-T1	1	Y	n7-18c22/23o138-166i
GPALN_002233-T1	0	Y	n4-11c21/22o
GPALN_013459-T1	0	Y	n5-16c22/23o
GPALN_006501-T1	0	Y	n12-25c30/31o
GPALN_010895-T2	0	Y	n4-15c20/21o
GPALN_007066-T1	0	Y	n6-13c17/18o
GPALN_004125-T1	1	Y	n4-13c18/19o73-96i
GPALN_006952-T1	0	Y	n15-25c30/31o
GPALN_009536-T1	0	Y	n3-13c17/18o
GPALN_008915-T1	0	Y	n2-13c18/19o
GPALN_009379-T1	1	Y	n9-19c24/25o232-262i
GPALN_002108-T1	0	Y	n7-18c24/25o
GPALN_001965-T1	1	Y	n2-13c18/19o302-323i
GPALN_001169-T1	0	Y	n2-13c18/19o
GPALN_008266-T1	0	Y	n2-12c17/18o
GPALN_004070-T1	0	Y	n2-12c19/20o
GPALN_014663-T1	9	Y	n10-20c25/26o35-53i65-83o103-122i143-164o205-224i236-257o263-285i320-340o352-375i
GPALN_015304-T1	0	Y	n4-15c20/21o
GPALN_011673-T1	0	Y	n7-15c23/24o
GPALN_005121-T1	0	Y	n17-28c36/37o
GPALN_007867-T1	0	Y	n3-14c18/19o
GPALN_006982-T1	0	Y	n2-9c13/14o
GPALN_009918-T1	0	Y	n2-20c24/25o
GPALN_010362-T1	0	Y	n2-12c17/18o
GPALN_007900-T1	1	Y	n7-18c23/24o460-482i
GPALN_002477-T1	0	Y	n3-13c18/19o
GPALN_005925-T1	0	Y	n2-9c16/17o
GPALN_013983-T1	2	Y	n3-14c18/19o1058-1077i1082-1102o
GPALN_012541-T1	0	Y	n24-35c43/44o
GPALN_004862-T1	0	Y	n3-15c20/21o
GPALN_014452-T1	0	Y	n4-11c16/17o
GPALN_013849-T1	1	Y	n6-17c22/23o982-1001i
GPALN_003369-T1	0	Y	n5-15c23/24o
GPALN_009137-T2	0	Y	n10-21c25/26o
GPALN_009670-T1	0	Y	n2-12c17/18o
GPALN_010794-T1	0	Y	n2-12c17/18o
GPALN_014846-T1	0	Y	n3-14c19/20o
GPALN_004557-T1	0	Y	n3-14c19/20o
GPALN_008784-T1	0	Y	n15-23c41/42o
GPALN_010427-T1	0	Y	n7-18c23/24o
GPALN_016202-T1	7	Y	n14-33c42/43o244-264i276-297o317-335i347-366o386-408i444-464o495-523i
GPALN_012415-T1	1	Y	n6-17c25/26o248-268i
GPALN_005101-T1	0	Y	n8-19c24/25o
GPALN_002147-T1	0	Y	n3-14c22/23o
GPALN_008357-T1	1	Y	n8-27c33/34o172-190i
GPALN_008484-T1	0	Y	n4-15c20/21o
GPALN_004075-T1	0	Y	n7-17c24/25o
GPALN_001111-T1	1	Y	n3-13c23/24o515-542i
GPALN_015910-T1	0	Y	n6-19c24/25o
GPALN_005973-T1	0	Y	n8-18c23/24o
GPALN_008603-T1	0	Y	n8-23c31/32o
GPALN_000589-T1	5	Y	n18-33c38/39o48-69i81-104o140-160i205-222o228-248i
GPALN_007851-T1	0	Y	n4-15c26/27o
GPALN_013283-T1	1	Y	n3-14c18/19o82-109i
GPALN_003409-T1	0	Y	n6-21c26/27o
GPALN_010824-T1	2	Y	n3-14c19/20o224-244i256-273o
GPALN_000164-T1	1	Y	n14-25c30/31o1371-1398i
GPALN_011305-T1	0	Y	n3-14c19/20o
GPALN_011710-T1	0	Y	n3-14c19/20o
GPALN_003560-T1	0	Y	n3-18c23/24o
GPALN_012417-T1	0	Y	n6-17c25/26o
GPALN_012209-T1	1	Y	n9-20c25/26o156-173i
GPALN_001564-T1	0	Y	n9-19c24/25o
GPALN_011458-T1	0	Y	n21-32c36/37o
GPALN_000177-T1	0	Y	n14-26c34/35o
GPALN_015782-T1	1	Y	n4-15c20/21o532-556i
GPALN_000475-T1	0	Y	n6-14c19/20o
GPALN_002125-T1	2	Y	n2-12c20/21o96-118i127-148o
GPALN_000182-T1	1	Y	n2-13c17/18o27-46i
GPALN_008086-T1	0	Y	n3-11c19/20o
GPALN_014099-T1	0	Y	n7-17c21/22o
GPALN_006854-T1	0	Y	n8-19c27/28o
GPALN_016180-T1	0	Y	n3-11c15/16o
GPALN_009943-T1	0	Y	n15-27c39/40o
GPALN_008016-T1	0	Y	n10-21c29/30o
GPALN_013400-T1	1	Y	n3-12c20/21o113-133i
GPALN_002184-T1	0	Y	n4-15c20/21o
GPALN_012969-T1	1	Y	n8-18c22/23o238-257i
GPALN_000992-T1	0	Y	n4-13c18/19o
GPALN_015927-T1	1	Y	n8-19c25/26o393-413i
GPALN_013544-T1	0	Y	n3-14c19/20o
GPALN_013502-T1	0	Y	n20-30c37/38o
GPALN_009580-T1	0	Y	n5-16c24/25o
GPALN_011775-T1	0	Y	n8-16c21/22o
GPALN_014294-T1	0	Y	n3-11c15/16o
GPALN_003705-T1	7	Y	n7-17c21/22o31-51i63-87o121-142i149-168o202-221i242-263o308-326i
GPALN_004414-T1	0	Y	n4-15c20/21o
GPALN_006804-T1	0	Y	n10-21c25/26o
GPALN_012296-T1	0	Y	n3-15c20/21o
GPALN_002781-T1	0	Y	n3-14c18/19o
GPALN_012211-T1	0	Y	n3-14c19/20o
GPALN_000509-T1	0	Y	n4-12c17/18o
GPALN_010126-T1	0	Y	n4-15c19/20o
GPALN_016351-T1	1	Y	n6-18c23/24o87-108i
GPALN_013601-T1	0	Y	n12-25c32/33o
GPALN_007964-T1	0	Y	n15-27c32/33o
GPALN_007584-T1	0	Y	n16-27c34/35o
GPALN_014672-T1	0	Y	n8-19c31/32o
GPALN_006239-T1	0	Y	n3-11c16/17o
GPALN_007083-T1	0	Y	n3-14c19/20o
GPALN_013348-T1	0	Y	n5-20c24/25o
GPALN_000601-T1	1	Y	n4-12c16/17o378-400i
GPALN_013442-T1	0	Y	n4-13c18/19o
GPALN_008126-T1	0	Y	n4-15c26/27o
GPALN_008666-T1	4	Y	n12-19c24/25o215-237i285-311o323-341i718-742o
GPALN_009055-T1	1	Y	n2-10c16/17o290-313i
GPALN_003352-T1	0	Y	n8-19c27/28o
GPALN_002420-T1	3	Y	n5-17c22/23o205-227i239-262o268-287i
GPALN_012839-T1	0	Y	n3-14c19/20o
GPALN_012822-T1	0	Y	n7-18c25/26o
GPALN_014747-T1	0	Y	n6-14c23/24o
GPALN_005140-T1	0	Y	n6-19c23/24o
GPALN_013701-T1	1	Y	n28-47c52/53o165-185i
GPALN_014050-T1	1	Y	n2-9c14/15o38-59i
GPALN_007181-T1	0	Y	n4-11c15/16o
GPALN_003986-T1	0	Y	n3-13c18/19o
GPALN_003478-T1	1	Y	n9-20c26/27o908-933i
GPALN_014697-T1	1	Y	n3-14c19/20o135-154i
GPALN_008574-T1	0	Y	n2-9c14/15o
GPALN_003006-T1	0	Y	n8-19c24/25o
GPALN_002883-T1	0	Y	n5-16c21/22o
GPALN_001167-T1	0	Y	n16-28c33/34o
GPALN_008370-T1	0	Y	n5-13c18/19o
GPALN_011495-T1	1	Y	n12-20c25/26o89-112i
GPALN_013598-T1	1	Y	n16-26c30/31o105-138i
GPALN_004560-T1	0	Y	n8-19c28/29o
GPALN_013634-T1	0	Y	n7-20c25/26o
GPALN_009552-T1	0	Y	n4-15c23/24o
GPALN_005072-T1	0	Y	n11-19c24/25o
GPALN_015875-T1	0	Y	n4-14c19/20o
GPALN_010534-T1	0	Y	n7-20c25/26o
GPALN_007028-T1	7	Y	n4-18c30/31o440-458i1040-1056o1076-1098i1105-1126o1146-1163i1170-1191o1197-1217i
GPALN_016065-T1	0	Y	n3-13c18/19o
GPALN_006507-T1	1	Y	n5-16c20/21o158-177i
GPALN_009640-T1	0	Y	n6-17c22/23o
GPALN_009953-T1	1	Y	n3-11c29/30o561-584i
GPALN_013844-T1	0	Y	n7-18c25/26o
GPALN_007035-T1	0	Y	n7-18c22/23o
GPALN_011856-T1	0	Y	n8-19c24/25o
GPALN_008085-T1	0	Y	n6-17c21/22o
GPALN_010645-T1	0	Y	n6-17c21/22o
GPALN_002990-T1	0	Y	n4-15c24/25o
GPALN_004881-T1	0	Y	n3-14c18/19o
GPALN_014268-T1	0	Y	n10-20c24/25o
GPALN_007529-T1	1	Y	n4-14c18/19o172-189i
GPALN_011184-T1	1	Y	n2-13c17/18o33-55i
GPALN_001069-T1	0	Y	n2-13c31/32o
GPALN_013233-T1	1	Y	n5-16c21/22o194-213i
GPALN_005728-T1	0	Y	n4-15c20/21o
GPALN_003679-T1	0	Y	n4-22c30/31o
GPALN_007884-T2	0	Y	n12-22c30/31o
GPALN_007325-T1	0	Y	n4-11c16/17o
GPALN_005057-T1	0	Y	n4-14c19/20o
GPALN_015791-T1	0	Y	n9-20c38/39o
GPALN_014061-T2	0	Y	n2-13c21/22o
GPALN_004572-T1	0	Y	n5-15c20/21o
GPALN_012034-T1	8	Y	n4-12c16/17o26-46i53-74o86-106i113-132o163-184i205-226o266-285i306-326o
GPALN_013066-T1	11	Y	n5-13c18/19o34-55i76-96o116-137i144-161o173-199i211-229o241-264i395-415o421-437i444-463o529-554i
GPALN_008151-T1	0	Y	n6-16c21/22o
GPALN_012410-T1	0	Y	n3-15c20/21o
GPALN_000972-T1	0	Y	n10-21c25/26o""".split("\n")

#############################################
##########################################
#########
signalP = """GPALN_000414-T1
GPALN_000424-T1
GPALN_000429-T1
GPALN_000430-T1
GPALN_000433-T1
GPALN_000434-T1
GPALN_000435-T1
GPALN_000438-T1
GPALN_000452-T1
GPALN_000453-T1
GPALN_000455-T1
GPALN_000484-T1
GPALN_000495-T1
GPALN_000509-T1
GPALN_000510-T1
GPALN_000517-T1
GPALN_000526-T1
GPALN_000527-T1
GPALN_000536-T1
GPALN_000538-T1
GPALN_000543-T1
GPALN_000555-T1
GPALN_000556-T1
GPALN_000558-T1
GPALN_000562-T1
GPALN_000566-T1
GPALN_000569-T1
GPALN_000576-T1
GPALN_000601-T1
GPALN_000608-T1
GPALN_000621-T1
GPALN_000629-T1
GPALN_000654-T1
GPALN_000680-T1
GPALN_000692-T1
GPALN_000698-T1
GPALN_000702-T1
GPALN_000702-T2
GPALN_000704-T1
GPALN_000705-T1
GPALN_000707-T1
GPALN_000767-T1
GPALN_000771-T1
GPALN_000778-T1
GPALN_000790-T1
# name         
GPALN_009357-T1
GPALN_009358-T1
GPALN_009364-T1
GPALN_009367-T1
GPALN_009374-T1
GPALN_009377-T1
GPALN_009377-T2
GPALN_009379-T1
GPALN_009380-T1
GPALN_009406-T1
GPALN_009412-T1
GPALN_009431-T1
GPALN_009441-T1
GPALN_009441-T2
GPALN_009441-T3
GPALN_009443-T1
GPALN_009444-T1
GPALN_009458-T1
GPALN_009466-T1
GPALN_009477-T1
GPALN_009490-T1
GPALN_009492-T1
GPALN_009493-T1
GPALN_009497-T1
GPALN_009498-T1
GPALN_009500-T1
GPALN_009502-T1
GPALN_009505-T1
GPALN_009509-T1
GPALN_009510-T1
GPALN_009512-T1
GPALN_009513-T1
GPALN_009532-T1
GPALN_009536-T1
GPALN_009537-T1
GPALN_009550-T1
GPALN_009552-T1
GPALN_009555-T2
GPALN_009556-T1
GPALN_009566-T1
GPALN_009568-T1
GPALN_009571-T1
GPALN_009571-T2
GPALN_009573-T1
GPALN_009579-T1
GPALN_009580-T1
GPALN_009582-T1
GPALN_009583-T1
GPALN_009584-T1
GPALN_009585-T1
GPALN_009586-T1
GPALN_009589-T1
GPALN_009590-T1
GPALN_009593-T1
GPALN_009609-T1
GPALN_009610-T1
GPALN_009613-T1
GPALN_009615-T1
GPALN_009618-T1
GPALN_009621-T1
GPALN_009622-T1
GPALN_009625-T1
GPALN_009627-T1
GPALN_009628-T1
GPALN_009629-T1
GPALN_009633-T1
GPALN_009636-T1
GPALN_009637-T1
GPALN_009638-T1
GPALN_009639-T1
GPALN_009640-T1
GPALN_009641-T1
GPALN_009648-T1
GPALN_009649-T1
GPALN_009650-T1
GPALN_009660-T1
GPALN_009663-T1
GPALN_009669-T1
GPALN_009670-T1
GPALN_009675-T1
GPALN_009682-T1
GPALN_009685-T1
GPALN_009695-T1
GPALN_009711-T1
GPALN_009716-T1
GPALN_009717-T1
GPALN_009721-T1
GPALN_009721-T2
GPALN_009723-T1
GPALN_009730-T1
GPALN_009731-T1
# name         
GPALN_001640-T1
GPALN_001641-T1
GPALN_001643-T1
GPALN_001668-T1
GPALN_001738-T1
GPALN_001748-T1
GPALN_001748-T2
GPALN_001751-T1
GPALN_001758-T1
GPALN_001759-T1
GPALN_001760-T1
GPALN_001780-T1
GPALN_001796-T1
GPALN_001796-T2
GPALN_001823-T1
GPALN_001828-T1
GPALN_001845-T2
GPALN_001846-T1
GPALN_001846-T2
GPALN_001847-T1
GPALN_001849-T1
GPALN_001862-T1
GPALN_001874-T1
GPALN_001881-T1
GPALN_001883-T1
GPALN_001887-T1
GPALN_001893-T1
GPALN_001899-T1
GPALN_001911-T1
GPALN_001912-T1
GPALN_001923-T1
GPALN_001924-T1
GPALN_001938-T1
GPALN_001951-T1
GPALN_001953-T1
GPALN_001954-T1
GPALN_001955-T1
GPALN_001965-T1
GPALN_001966-T1
GPALN_001972-T1
GPALN_001982-T1
GPALN_001989-T1
GPALN_002001-T1
GPALN_002018-T2
GPALN_002028-T1
# name         
GPALN_015155-T1
GPALN_015156-T1
GPALN_015157-T1
GPALN_015161-T1
GPALN_015163-T1
GPALN_015169-T1
GPALN_015172-T1
GPALN_015174-T1
GPALN_015174-T2
GPALN_015175-T1
GPALN_015177-T1
GPALN_015178-T1
GPALN_015179-T1
GPALN_015180-T1
GPALN_015181-T1
GPALN_015182-T1
GPALN_015183-T1
GPALN_015184-T1
GPALN_015186-T1
GPALN_015188-T1
GPALN_015193-T1
GPALN_015209-T1
GPALN_015211-T1
GPALN_015217-T1
GPALN_015218-T1
GPALN_015222-T1
GPALN_015223-T1
GPALN_015224-T1
GPALN_015227-T1
GPALN_015230-T1
GPALN_015232-T1
GPALN_015236-T1
GPALN_015237-T1
GPALN_015243-T1
GPALN_015244-T1
GPALN_015245-T1
GPALN_015246-T1
GPALN_015248-T1
GPALN_015251-T1
GPALN_015252-T1
GPALN_015253-T1
GPALN_015254-T1
GPALN_015258-T1
GPALN_015260-T1
GPALN_015267-T1
GPALN_015272-T1
GPALN_015279-T1
GPALN_015280-T1
GPALN_015281-T1
GPALN_015285-T1
GPALN_015286-T1
GPALN_015287-T1
GPALN_015288-T1
GPALN_015291-T1
GPALN_015292-T1
GPALN_015294-T1
GPALN_015295-T1
GPALN_015296-T1
GPALN_015297-T1
GPALN_015298-T1
GPALN_015299-T1
GPALN_015300-T1
GPALN_015301-T1
GPALN_015302-T1
GPALN_015304-T1
GPALN_015307-T1
GPALN_015309-T1
GPALN_015314-T1
GPALN_015318-T1
GPALN_015319-T1
GPALN_015324-T1
GPALN_015326-T1
GPALN_015332-T1
GPALN_015334-T1
GPALN_015347-T1
GPALN_015354-T1
GPALN_015355-T1
GPALN_015357-T1
GPALN_015358-T1
GPALN_015359-T1
GPALN_015361-T1
GPALN_015364-T1
GPALN_015365-T1
GPALN_015367-T1
GPALN_015372-T1
GPALN_015372-T2
GPALN_015381-T1
GPALN_015408-T1
GPALN_015415-T1
GPALN_015418-T1
GPALN_015423-T1
GPALN_015425-T1
GPALN_015436-T1
GPALN_015446-T1
GPALN_015447-T1
GPALN_015450-T1
GPALN_015469-T1
GPALN_015495-T1
GPALN_015499-T1
GPALN_015499-T2
GPALN_015505-T1
GPALN_015531-T1
GPALN_015532-T1
GPALN_015542-T1
# name         
GPALN_000018-T1
GPALN_000040-T1
GPALN_000045-T1
GPALN_000048-T1
GPALN_000071-T1
GPALN_000079-T1
GPALN_000109-T1
GPALN_000112-T1
GPALN_000115-T1
GPALN_000123-T1
GPALN_000124-T1
GPALN_000149-T1
GPALN_000162-T1
GPALN_000163-T1
GPALN_000164-T1
GPALN_000167-T1
GPALN_000171-T1
GPALN_000176-T1
GPALN_000181-T1
GPALN_000196-T1
GPALN_000253-T1
GPALN_000262-T1
GPALN_000274-T1
GPALN_000311-T1
GPALN_000323-T1
GPALN_000341-T1
GPALN_000342-T1
GPALN_000343-T1
GPALN_000360-T1
GPALN_000361-T1
GPALN_000371-T1
GPALN_000375-T1
GPALN_000377-T1
GPALN_000383-T1
GPALN_000385-T1
GPALN_000393-T1
GPALN_000394-T1
GPALN_000395-T1
GPALN_000397-T1
# name         
GPALN_011856-T1
GPALN_011857-T1
GPALN_011858-T1
GPALN_011865-T1
GPALN_011865-T2
GPALN_011867-T1
GPALN_011875-T1
GPALN_011879-T1
GPALN_011881-T1
GPALN_011884-T1
GPALN_011891-T1
GPALN_011891-T2
GPALN_011894-T1
GPALN_011917-T1
GPALN_011921-T1
GPALN_011927-T1
GPALN_011955-T1
GPALN_011956-T1
GPALN_011957-T1
GPALN_011972-T1
GPALN_011989-T1
GPALN_012004-T1
GPALN_012005-T1
GPALN_012007-T1
GPALN_012010-T1
GPALN_012021-T1
GPALN_012022-T1
GPALN_012025-T1
GPALN_012056-T1
GPALN_012062-T1
GPALN_012064-T1
GPALN_012067-T1
GPALN_012077-T1
GPALN_012096-T1
GPALN_012099-T1
GPALN_012122-T1
GPALN_012123-T1
GPALN_012124-T1
GPALN_012127-T1
GPALN_012134-T1
GPALN_012140-T1
GPALN_012150-T1
GPALN_012183-T1
GPALN_012186-T1
GPALN_012189-T1
GPALN_012205-T1
GPALN_012207-T1
GPALN_012209-T1
GPALN_012215-T1
GPALN_012232-T1
GPALN_012234-T1
# name         
GPALN_004064-T1
GPALN_004065-T1
GPALN_004070-T1
GPALN_004071-T1
GPALN_004073-T1
GPALN_004075-T1
GPALN_004076-T1
GPALN_004078-T1
GPALN_004105-T1
GPALN_004107-T1
GPALN_004117-T1
GPALN_004117-T2
GPALN_004125-T1
GPALN_004130-T1
GPALN_004131-T1
GPALN_004166-T1
GPALN_004167-T1
GPALN_004168-T1
GPALN_004174-T1
GPALN_004178-T1
GPALN_004179-T1
GPALN_004180-T1
GPALN_004194-T1
GPALN_004195-T1
GPALN_004197-T1
GPALN_004203-T1
GPALN_004209-T1
GPALN_004239-T1
GPALN_004239-T2
GPALN_004253-T1
GPALN_004265-T1
GPALN_004288-T1
GPALN_004292-T1
GPALN_004293-T1
GPALN_004294-T1
GPALN_004306-T1
GPALN_004308-T1
GPALN_004316-T1
GPALN_004342-T1
GPALN_004346-T1
GPALN_004347-T1
GPALN_004349-T1
GPALN_004350-T1
GPALN_004355-T1
GPALN_004368-T1
GPALN_004369-T1
GPALN_004373-T1
GPALN_004377-T1
GPALN_004380-T1
GPALN_004382-T1
GPALN_004384-T1
GPALN_004398-T1
GPALN_004399-T1
GPALN_004410-T1
GPALN_004411-T1
GPALN_004414-T1
GPALN_004418-T1
GPALN_004425-T1
GPALN_004440-T1
GPALN_004447-T1
GPALN_004450-T1
GPALN_004458-T1
GPALN_004460-T1
GPALN_004463-T1
GPALN_004464-T1
GPALN_004470-T1
GPALN_004479-T1
GPALN_004480-T1
GPALN_004481-T1
# name         
GPALN_014730-T1
GPALN_014738-T1
GPALN_014745-T1
GPALN_014747-T1
GPALN_014750-T1
GPALN_014753-T1
GPALN_014755-T1
GPALN_014759-T1
GPALN_014775-T1
GPALN_014785-T1
GPALN_014793-T1
GPALN_014799-T1
GPALN_014800-T1
GPALN_014802-T1
GPALN_014831-T1
GPALN_014835-T1
GPALN_014842-T1
GPALN_014843-T1
GPALN_014845-T1
GPALN_014846-T1
GPALN_014851-T1
GPALN_014857-T1
GPALN_014862-T1
GPALN_014865-T1
GPALN_014866-T1
GPALN_014867-T1
GPALN_014868-T1
GPALN_014870-T1
GPALN_014874-T1
GPALN_014876-T1
GPALN_014879-T1
GPALN_014881-T1
GPALN_014883-T1
GPALN_014884-T1
GPALN_014885-T1
GPALN_014890-T1
GPALN_014893-T1
GPALN_014897-T1
GPALN_014903-T1
GPALN_014904-T1
GPALN_014919-T1
GPALN_014924-T1
GPALN_014929-T1
GPALN_014935-T1
GPALN_014938-T1
GPALN_014939-T1
GPALN_014951-T1
GPALN_014952-T1
GPALN_014957-T1
GPALN_014960-T1
GPALN_014962-T1
GPALN_014964-T1
GPALN_014966-T1
GPALN_014967-T1
GPALN_014971-T1
GPALN_014975-T1
GPALN_014976-T1
GPALN_015006-T1
GPALN_015012-T1
GPALN_015013-T1
GPALN_015025-T1
GPALN_015037-T1
GPALN_015046-T1
GPALN_015046-T2
GPALN_015046-T3
GPALN_015061-T1
GPALN_015073-T1
GPALN_015075-T1
GPALN_015094-T1
GPALN_015100-T1
GPALN_015116-T1
GPALN_015121-T1
GPALN_015123-T1
# name         
GPALN_012664-T1
GPALN_012667-T1
GPALN_012702-T1
GPALN_012703-T1
GPALN_012709-T1
GPALN_012718-T1
GPALN_012727-T1
GPALN_012731-T1
GPALN_012734-T1
GPALN_012760-T1
GPALN_012762-T1
GPALN_012766-T1
GPALN_012791-T1
GPALN_012815-T1
GPALN_012822-T1
GPALN_012828-T1
GPALN_012829-T1
GPALN_012839-T1
GPALN_012841-T1
GPALN_012842-T1
GPALN_012843-T1
GPALN_012846-T1
GPALN_012854-T1
GPALN_012873-T1
GPALN_012887-T1
GPALN_012894-T1
GPALN_012897-T1
GPALN_012912-T1
GPALN_012921-T1
GPALN_012925-T1
GPALN_012947-T1
GPALN_012955-T1
GPALN_012956-T1
GPALN_012964-T1
GPALN_012969-T1
GPALN_012973-T1
GPALN_012979-T1
GPALN_012981-T1
GPALN_012989-T1
GPALN_012998-T1
GPALN_012998-T2
GPALN_013004-T1
GPALN_013005-T1
GPALN_013039-T1
# name         
GPALN_013883-T1
GPALN_013913-T1
GPALN_013923-T1
GPALN_013937-T1
GPALN_013939-T1
GPALN_013941-T1
GPALN_013983-T1
GPALN_013990-T1
GPALN_014001-T1
GPALN_014002-T1
GPALN_014005-T1
GPALN_014006-T1
GPALN_014019-T1
GPALN_014028-T1
GPALN_014032-T1
GPALN_014033-T1
GPALN_014034-T1
GPALN_014038-T1
GPALN_014042-T1
GPALN_014044-T1
GPALN_014045-T1
GPALN_014052-T1
GPALN_014060-T1
GPALN_014061-T1
GPALN_014061-T2
GPALN_014075-T1
GPALN_014079-T1
GPALN_014080-T1
GPALN_014086-T1
GPALN_014087-T1
GPALN_014091-T1
GPALN_014094-T1
GPALN_014095-T1
GPALN_014099-T1
GPALN_014102-T1
GPALN_014107-T1
GPALN_014121-T1
GPALN_014122-T1
GPALN_014126-T1
GPALN_014133-T1
GPALN_014145-T1
GPALN_014146-T1
GPALN_014147-T1
GPALN_014150-T1
GPALN_014151-T1
GPALN_014172-T1
GPALN_014182-T1
GPALN_014189-T1
GPALN_014191-T1
GPALN_014226-T1
GPALN_014231-T1
GPALN_014232-T1
GPALN_014235-T1
GPALN_014237-T1
GPALN_014238-T1
GPALN_014261-T1
GPALN_014264-T1
GPALN_014268-T1
GPALN_014271-T1
GPALN_014279-T1
# name         
GPALN_006126-T1
GPALN_006154-T1
GPALN_006162-T1
GPALN_006176-T1
GPALN_006201-T1
GPALN_006223-T1
GPALN_006252-T1
GPALN_006261-T1
GPALN_006283-T1
GPALN_006286-T1
GPALN_006304-T1
GPALN_006313-T1
GPALN_006318-T1
GPALN_006319-T1
GPALN_006328-T1
GPALN_006343-T1
GPALN_006346-T1
GPALN_006369-T1
GPALN_006381-T1
GPALN_006406-T1
GPALN_006413-T1
GPALN_006425-T1
GPALN_006433-T1
GPALN_006440-T1
GPALN_006445-T1
GPALN_006471-T1
GPALN_006482-T1
GPALN_006486-T1
GPALN_006490-T1
GPALN_006501-T1
GPALN_006507-T1
# name         
GPALN_013480-T1
GPALN_013481-T1
GPALN_013482-T1
GPALN_013490-T1
GPALN_013496-T1
GPALN_013504-T1
GPALN_013504-T2
GPALN_013508-T1
GPALN_013514-T1
GPALN_013515-T1
GPALN_013528-T1
GPALN_013533-T1
GPALN_013535-T1
GPALN_013542-T1
GPALN_013544-T1
GPALN_013550-T1
GPALN_013561-T1
GPALN_013563-T1
GPALN_013575-T1
GPALN_013577-T1
GPALN_013578-T1
GPALN_013580-T1
GPALN_013581-T1
GPALN_013598-T1
GPALN_013613-T1
GPALN_013623-T1
GPALN_013634-T1
GPALN_013639-T1
GPALN_013660-T1
GPALN_013662-T1
GPALN_013702-T1
GPALN_013702-T2
GPALN_013707-T1
GPALN_013720-T1
GPALN_013737-T1
GPALN_013749-T1
GPALN_013767-T1
GPALN_013771-T1
GPALN_013781-T1
GPALN_013782-T1
GPALN_013782-T2
GPALN_013784-T1
GPALN_013793-T1
GPALN_013816-T1
GPALN_013819-T1
GPALN_013821-T1
GPALN_013828-T1
GPALN_013829-T1
GPALN_013844-T1
GPALN_013849-T1
GPALN_013858-T1
# name         
GPALN_000842-T1
GPALN_000858-T1
GPALN_000858-T2
GPALN_000870-T1
GPALN_000889-T1
GPALN_000894-T1
GPALN_000922-T1
GPALN_000924-T1
GPALN_000928-T1
GPALN_000937-T1
GPALN_000945-T1
GPALN_000948-T1
GPALN_000949-T1
GPALN_000953-T1
GPALN_000972-T1
GPALN_000982-T1
GPALN_000992-T1
GPALN_000994-T1
GPALN_001002-T1
GPALN_001004-T1
GPALN_001008-T1
GPALN_001012-T1
GPALN_001013-T1
GPALN_001014-T1
GPALN_001017-T1
GPALN_001018-T1
GPALN_001018-T2
GPALN_001026-T1
GPALN_001041-T1
GPALN_001042-T1
GPALN_001055-T1
GPALN_001062-T1
GPALN_001063-T1
GPALN_001069-T1
GPALN_001091-T1
GPALN_001093-T1
GPALN_001094-T1
GPALN_001101-T1
GPALN_001106-T1
GPALN_001110-T1
GPALN_001111-T1
GPALN_001115-T1
GPALN_001121-T1
GPALN_001134-T1
GPALN_001139-T1
GPALN_001144-T1
GPALN_001145-T1
GPALN_001146-T1
GPALN_001147-T1
GPALN_001149-T1
GPALN_001150-T1
GPALN_001151-T1
GPALN_001153-T1
GPALN_001161-T1
GPALN_001162-T1
GPALN_001169-T1
GPALN_001174-T1
GPALN_001179-T1
GPALN_001182-T1
GPALN_001186-T1
GPALN_001188-T1
GPALN_001190-T1
GPALN_001202-T1
GPALN_001203-T1
GPALN_001207-T1
GPALN_001208-T1
# name         
GPALN_008553-T1
GPALN_008558-T1
GPALN_008584-T1
GPALN_008597-T1
GPALN_008601-T1
GPALN_008603-T1
GPALN_008619-T1
GPALN_008620-T1
GPALN_008633-T1
GPALN_008643-T1
GPALN_008658-T1
GPALN_008664-T1
GPALN_008666-T1
GPALN_008673-T1
GPALN_008680-T1
GPALN_008681-T1
GPALN_008681-T2
GPALN_008695-T1
GPALN_008701-T1
GPALN_008738-T1
GPALN_008744-T1
GPALN_008745-T1
GPALN_008762-T1
GPALN_008764-T1
GPALN_008771-T1
GPALN_008791-T1
GPALN_008808-T1
GPALN_008814-T1
GPALN_008815-T1
GPALN_008816-T1
GPALN_008824-T1
GPALN_008839-T1
GPALN_008852-T1
GPALN_008859-T1
GPALN_008860-T1
GPALN_008866-T1
GPALN_008867-T1
GPALN_008878-T1
GPALN_008881-T1
GPALN_008886-T1
GPALN_008908-T1
GPALN_008910-T1
GPALN_008914-T1
GPALN_008915-T1
# name         
GPALN_002450-T1
GPALN_002453-T1
GPALN_002455-T1
GPALN_002460-T1
GPALN_002464-T1
GPALN_002466-T1
GPALN_002471-T1
GPALN_002472-T1
GPALN_002477-T1
GPALN_002478-T1
GPALN_002481-T1
GPALN_002482-T1
GPALN_002492-T1
GPALN_002493-T1
GPALN_002494-T1
GPALN_002500-T1
GPALN_002508-T1
GPALN_002509-T1
GPALN_002515-T1
GPALN_002516-T1
GPALN_002519-T1
GPALN_002521-T1
GPALN_002527-T1
GPALN_002542-T1
GPALN_002546-T1
GPALN_002547-T1
GPALN_002550-T1
GPALN_002552-T1
GPALN_002554-T1
GPALN_002562-T1
GPALN_002574-T1
GPALN_002578-T1
GPALN_002593-T1
GPALN_002595-T1
GPALN_002600-T1
GPALN_002601-T1
GPALN_002605-T1
GPALN_002641-T1
GPALN_002651-T1
GPALN_002656-T1
GPALN_002666-T1
GPALN_002671-T1
GPALN_002674-T1
GPALN_002678-T1
GPALN_002681-T1
GPALN_002682-T1
GPALN_002687-T1
GPALN_002690-T1
GPALN_002692-T1
GPALN_002695-T1
GPALN_002707-T1
GPALN_002710-T1
GPALN_002713-T1
GPALN_002729-T1
GPALN_002750-T1
GPALN_002759-T1
GPALN_002760-T1
GPALN_002762-T1
GPALN_002769-T1
GPALN_002770-T1
GPALN_002775-T1
GPALN_002778-T1
GPALN_002779-T1
GPALN_002781-T1
GPALN_002784-T1
GPALN_002788-T1
GPALN_002791-T1
GPALN_002802-T1
# name         
GPALN_003663-T1
GPALN_003679-T1
GPALN_003683-T1
GPALN_003693-T1
GPALN_003712-T1
GPALN_003714-T1
GPALN_003719-T1
GPALN_003741-T1
GPALN_003743-T1
GPALN_003757-T1
GPALN_003757-T2
GPALN_003759-T1
GPALN_003760-T1
GPALN_003769-T1
GPALN_003771-T1
GPALN_003772-T1
GPALN_003786-T1
GPALN_003792-T1
GPALN_003793-T1
GPALN_003794-T1
GPALN_003795-T1
GPALN_003797-T1
GPALN_003803-T1
GPALN_003809-T1
GPALN_003816-T1
GPALN_003819-T1
GPALN_003820-T1
GPALN_003824-T1
GPALN_003825-T1
GPALN_003828-T1
GPALN_003829-T1
GPALN_003831-T1
GPALN_003832-T1
GPALN_003835-T1
GPALN_003836-T1
GPALN_003837-T1
GPALN_003838-T1
GPALN_003839-T1
GPALN_003844-T1
GPALN_003846-T1
GPALN_003850-T1
GPALN_003851-T1
GPALN_003852-T1
GPALN_003860-T1
GPALN_003867-T1
GPALN_003873-T1
GPALN_003876-T1
GPALN_003882-T1
GPALN_003883-T1
GPALN_003890-T1
GPALN_003891-T1
GPALN_003895-T1
GPALN_003897-T1
GPALN_003903-T1
GPALN_003905-T1
GPALN_003906-T1
GPALN_003907-T2
GPALN_003908-T1
GPALN_003910-T1
GPALN_003911-T1
GPALN_003912-T1
GPALN_003913-T1
GPALN_003917-T1
GPALN_003918-T1
GPALN_003925-T1
GPALN_003937-T1
GPALN_003941-T1
GPALN_003942-T1
GPALN_003946-T1
GPALN_003949-T1
GPALN_003950-T1
GPALN_003951-T1
GPALN_003952-T1
GPALN_003953-T1
GPALN_003954-T1
GPALN_003955-T1
GPALN_003956-T1
GPALN_003957-T1
GPALN_003959-T1
GPALN_003961-T1
GPALN_003969-T1
GPALN_003970-T1
GPALN_003971-T1
GPALN_003974-T1
GPALN_003975-T1
GPALN_003977-T1
GPALN_003980-T1
GPALN_003983-T1
GPALN_003986-T1
GPALN_003990-T1
GPALN_003997-T1
GPALN_003999-T1
GPALN_004000-T1
GPALN_004007-T1
GPALN_004008-T1
GPALN_004009-T1
GPALN_004011-T1
GPALN_004013-T1
GPALN_004014-T1
GPALN_004017-T1
GPALN_004018-T1
GPALN_004019-T1
GPALN_004019-T2
GPALN_004020-T1
GPALN_004023-T1
GPALN_004025-T1
GPALN_004026-T1
GPALN_004031-T1
GPALN_004045-T1
GPALN_004050-T1
GPALN_004056-T1
GPALN_004057-T1
GPALN_004059-T1
# name         
GPALN_005309-T1
GPALN_005326-T1
GPALN_005373-T1
GPALN_005393-T1
GPALN_005399-T1
GPALN_005401-T1
GPALN_005416-T1
GPALN_005418-T1
GPALN_005419-T1
GPALN_005439-T1
GPALN_005443-T1
GPALN_005451-T1
GPALN_005452-T1
GPALN_005489-T1
GPALN_005501-T1
GPALN_005531-T1
GPALN_005537-T1
GPALN_005538-T1
GPALN_005553-T1
GPALN_005554-T1
GPALN_005574-T1
GPALN_005577-T1
GPALN_005591-T1
GPALN_005598-T1
GPALN_005599-T1
GPALN_005609-T1
GPALN_005611-T1
GPALN_005612-T1
GPALN_005613-T1
GPALN_005617-T1
GPALN_005637-T1
GPALN_005645-T1
GPALN_005679-T1
# name         
GPALN_008978-T1
GPALN_008999-T1
GPALN_009008-T1
GPALN_009015-T1
GPALN_009029-T1
GPALN_009029-T2
GPALN_009038-T1
GPALN_009039-T1
GPALN_009055-T1
GPALN_009056-T1
GPALN_009074-T1
GPALN_009084-T1
GPALN_009091-T1
GPALN_009103-T1
GPALN_009111-T1
GPALN_009147-T1
GPALN_009156-T1
GPALN_009205-T1
GPALN_009225-T1
GPALN_009233-T1
GPALN_009253-T1
GPALN_009273-T1
GPALN_009281-T1
GPALN_009291-T1
GPALN_009292-T1
GPALN_009302-T1
GPALN_009314-T1
GPALN_009323-T1
GPALN_009325-T1
GPALN_009333-T4
GPALN_009335-T1
# name         
GPALN_007748-T1
GPALN_007756-T1
GPALN_007795-T1
GPALN_007796-T1
GPALN_007805-T1
GPALN_007806-T1
GPALN_007811-T1
GPALN_007820-T1
GPALN_007837-T1
GPALN_007838-T1
GPALN_007844-T1
GPALN_007846-T1
GPALN_007848-T1
GPALN_007850-T1
GPALN_007851-T1
GPALN_007859-T1
GPALN_007864-T1
GPALN_007866-T1
GPALN_007867-T1
GPALN_007871-T1
GPALN_007884-T1
GPALN_007884-T2
GPALN_007886-T1
GPALN_007899-T1
GPALN_007900-T1
GPALN_007905-T1
GPALN_007916-T1
GPALN_007919-T1
GPALN_007921-T1
GPALN_007923-T1
GPALN_007934-T1
GPALN_007962-T1
GPALN_007964-T1
GPALN_007974-T1
GPALN_007976-T1
GPALN_007981-T1
GPALN_007985-T1
GPALN_007990-T1
GPALN_007995-T1
GPALN_008002-T1
GPALN_008006-T1
GPALN_008016-T1
GPALN_008042-T1
GPALN_008074-T1
GPALN_008078-T1
GPALN_008079-T1
GPALN_008081-T1
GPALN_008082-T1
GPALN_008083-T1
GPALN_008085-T1
GPALN_008086-T1
GPALN_008086-T2
GPALN_008088-T1
GPALN_008089-T1
GPALN_008097-T1
GPALN_008098-T1
GPALN_008100-T1
GPALN_008101-T1
GPALN_008102-T1
GPALN_008104-T1
GPALN_008105-T1
GPALN_008106-T1
GPALN_008107-T1
GPALN_008107-T2
GPALN_008108-T1
GPALN_008110-T1
GPALN_008111-T1
GPALN_008112-T1
GPALN_008113-T1
GPALN_008126-T1
GPALN_008130-T1
GPALN_008135-T1
GPALN_008137-T1
GPALN_008141-T1
GPALN_008144-T1
GPALN_008145-T1
# name         
GPALN_003264-T1
GPALN_003283-T1
GPALN_003306-T1
GPALN_003308-T1
GPALN_003309-T1
GPALN_003326-T1
GPALN_003328-T1
GPALN_003329-T1
GPALN_003340-T1
GPALN_003352-T1
GPALN_003357-T1
GPALN_003368-T1
GPALN_003369-T1
GPALN_003370-T1
GPALN_003371-T1
GPALN_003375-T1
GPALN_003376-T1
GPALN_003377-T1
GPALN_003379-T1
GPALN_003380-T1
GPALN_003381-T1
GPALN_003382-T1
GPALN_003383-T1
GPALN_003396-T1
GPALN_003398-T1
GPALN_003402-T1
GPALN_003405-T1
GPALN_003409-T1
GPALN_003415-T1
GPALN_003416-T1
GPALN_003416-T2
GPALN_003420-T1
GPALN_003422-T1
GPALN_003426-T1
GPALN_003427-T1
GPALN_003430-T1
GPALN_003431-T1
GPALN_003438-T1
GPALN_003445-T1
GPALN_003447-T1
GPALN_003449-T1
GPALN_003454-T1
GPALN_003462-T1
GPALN_003463-T1
GPALN_003468-T1
GPALN_003470-T1
GPALN_003478-T1
GPALN_003481-T1
GPALN_003486-T1
GPALN_003488-T1
GPALN_003547-T1
GPALN_003557-T1
GPALN_003566-T1
GPALN_003577-T1
GPALN_003580-T1
GPALN_003588-T1
GPALN_003596-T1
GPALN_003606-T1
GPALN_003621-T1
GPALN_003632-T1
GPALN_003641-T1
GPALN_003643-T1
GPALN_003649-T1
# name         
GPALN_005697-T1
GPALN_005711-T1
GPALN_005712-T1
GPALN_005717-T1
GPALN_005725-T1
GPALN_005728-T1
GPALN_005731-T1
GPALN_005732-T1
GPALN_005734-T1
GPALN_005736-T1
GPALN_005738-T1
GPALN_005741-T1
GPALN_005742-T1
GPALN_005743-T1
GPALN_005745-T1
GPALN_005749-T1
GPALN_005750-T1
GPALN_005751-T1
GPALN_005752-T1
GPALN_005753-T1
GPALN_005754-T1
GPALN_005758-T1
GPALN_005765-T1
GPALN_005766-T1
GPALN_005769-T1
GPALN_005773-T1
GPALN_005776-T1
GPALN_005777-T1
GPALN_005789-T1
GPALN_005799-T1
GPALN_005801-T1
GPALN_005802-T1
GPALN_005809-T1
GPALN_005812-T1
GPALN_005815-T1
GPALN_005818-T1
GPALN_005820-T1
GPALN_005821-T1
GPALN_005823-T1
GPALN_005829-T1
GPALN_005846-T1
GPALN_005866-T1
GPALN_005867-T1
GPALN_005868-T1
GPALN_005870-T1
GPALN_005872-T1
GPALN_005882-T1
GPALN_005889-T1
GPALN_005893-T1
GPALN_005896-T1
GPALN_005897-T1
GPALN_005901-T1
GPALN_005902-T1
GPALN_005903-T1
GPALN_005904-T1
GPALN_005909-T1
GPALN_005910-T1
GPALN_005911-T1
GPALN_005912-T1
GPALN_005915-T1
GPALN_005916-T1
GPALN_005918-T1
GPALN_005924-T1
GPALN_005925-T1
GPALN_005926-T1
GPALN_005927-T1
GPALN_005928-T1
GPALN_005930-T1
GPALN_005939-T1
GPALN_005953-T1
GPALN_005969-T1
GPALN_005970-T1
GPALN_005973-T1
GPALN_006002-T1
GPALN_006005-T1
GPALN_006022-T1
GPALN_006023-T1
GPALN_006026-T1
GPALN_006027-T1
GPALN_006028-T1
GPALN_006029-T1
GPALN_006031-T1
GPALN_006035-T1
GPALN_006038-T1
GPALN_006039-T2
GPALN_006041-T1
GPALN_006042-T1
GPALN_006043-T1
GPALN_006046-T1
GPALN_006047-T1
GPALN_006053-T1
GPALN_006055-T1
GPALN_006057-T1
GPALN_006059-T1
GPALN_006061-T1
GPALN_006063-T1
GPALN_006067-T1
GPALN_006072-T1
GPALN_006085-T1
GPALN_006086-T1
GPALN_006088-T1
GPALN_006088-T2
GPALN_006100-T1
GPALN_006102-T1
GPALN_006112-T1
# name         
GPALN_007357-T1
GPALN_007361-T1
GPALN_007362-T1
GPALN_007365-T1
GPALN_007368-T1
GPALN_007369-T1
GPALN_007371-T1
GPALN_007372-T1
GPALN_007374-T1
GPALN_007376-T1
GPALN_007380-T1
GPALN_007385-T1
GPALN_007386-T1
GPALN_007404-T1
GPALN_007409-T1
GPALN_007410-T1
GPALN_007411-T1
GPALN_007416-T1
GPALN_007422-T1
GPALN_007429-T1
GPALN_007433-T1
GPALN_007436-T1
GPALN_007443-T1
GPALN_007445-T1
GPALN_007468-T1
GPALN_007473-T1
GPALN_007473-T2
GPALN_007488-T1
GPALN_007498-T1
GPALN_007502-T1
GPALN_007510-T1
GPALN_007512-T1
GPALN_007514-T1
GPALN_007525-T1
GPALN_007529-T1
GPALN_007532-T1
GPALN_007536-T1
GPALN_007537-T1
GPALN_007539-T1
GPALN_007541-T1
GPALN_007543-T1
GPALN_007544-T1
GPALN_007548-T1
GPALN_007549-T1
GPALN_007550-T1
GPALN_007552-T1
GPALN_007554-T1
GPALN_007556-T1
GPALN_007557-T1
GPALN_007560-T1
GPALN_007561-T1
GPALN_007564-T1
GPALN_007565-T1
GPALN_007568-T1
GPALN_007574-T1
GPALN_007574-T2
GPALN_007580-T1
GPALN_007584-T1
GPALN_007586-T1
GPALN_007589-T1
GPALN_007592-T1
GPALN_007592-T2
GPALN_007592-T3
GPALN_007597-T1
GPALN_007603-T1
GPALN_007604-T1
GPALN_007605-T1
GPALN_007606-T1
GPALN_007616-T1
GPALN_007617-T1
GPALN_007619-T1
GPALN_007620-T1
GPALN_007622-T1
GPALN_007627-T1
GPALN_007633-T1
GPALN_007638-T1
GPALN_007643-T1
GPALN_007646-T1
GPALN_007647-T1
GPALN_007648-T1
GPALN_007670-T1
GPALN_007672-T1
GPALN_007673-T1
GPALN_007680-T1
GPALN_007681-T1
GPALN_007682-T1
GPALN_007683-T1
GPALN_007684-T1
GPALN_007685-T1
GPALN_007692-T1
GPALN_007693-T1
GPALN_007696-T1
GPALN_007697-T1
GPALN_007699-T1
GPALN_007704-T1
GPALN_007705-T1
GPALN_007706-T1
GPALN_007708-T1
GPALN_007710-T1
GPALN_007711-T1
GPALN_007715-T1
GPALN_007719-T1
GPALN_007720-T1
GPALN_007729-T1
GPALN_007731-T1
GPALN_007732-T1
GPALN_007737-T1
GPALN_007739-T1
GPALN_007740-T1
# name         
GPALN_008150-T1
GPALN_008151-T1
GPALN_008152-T1
GPALN_008154-T1
GPALN_008161-T1
GPALN_008176-T1
GPALN_008215-T1
GPALN_008228-T1
GPALN_008232-T1
GPALN_008245-T1
GPALN_008251-T1
GPALN_008266-T1
GPALN_008303-T1
GPALN_008308-T1
GPALN_008317-T1
GPALN_008323-T1
GPALN_008335-T1
GPALN_008335-T2
GPALN_008336-T1
GPALN_008339-T1
GPALN_008342-T1
GPALN_008356-T1
GPALN_008357-T1
GPALN_008366-T1
GPALN_008377-T1
GPALN_008379-T1
GPALN_008391-T1
GPALN_008420-T1
GPALN_008457-T1
GPALN_008474-T1
GPALN_008479-T1
GPALN_008484-T1
GPALN_008500-T1
GPALN_008501-T1
GPALN_008502-T1
GPALN_008510-T1
GPALN_008511-T1
GPALN_008515-T1
GPALN_008535-T1
GPALN_008542-T1
GPALN_008543-T1
# name         
GPALN_014292-T1
GPALN_014294-T1
GPALN_014304-T1
GPALN_014304-T2
GPALN_014308-T1
GPALN_014309-T1
GPALN_014324-T1
GPALN_014327-T1
GPALN_014329-T1
GPALN_014334-T1
GPALN_014351-T1
GPALN_014352-T1
GPALN_014354-T1
GPALN_014355-T1
GPALN_014357-T1
GPALN_014368-T1
GPALN_014369-T1
GPALN_014370-T1
GPALN_014371-T1
GPALN_014372-T1
GPALN_014376-T1
GPALN_014377-T1
GPALN_014378-T1
GPALN_014379-T1
GPALN_014381-T1
GPALN_014392-T1
GPALN_014395-T1
GPALN_014397-T1
GPALN_014398-T1
GPALN_014403-T1
GPALN_014408-T1
GPALN_014415-T1
GPALN_014416-T1
GPALN_014438-T1
GPALN_014459-T1
GPALN_014471-T1
GPALN_014477-T1
GPALN_014479-T1
GPALN_014498-T1
GPALN_014500-T1
GPALN_014501-T1
GPALN_014503-T1
GPALN_014523-T1
GPALN_014529-T1
GPALN_014530-T1
GPALN_014539-T1
GPALN_014549-T1
GPALN_014552-T1
GPALN_014567-T1
GPALN_014568-T1
GPALN_014569-T1
GPALN_014571-T1
GPALN_014575-T1
GPALN_014576-T1
GPALN_014583-T1
GPALN_014639-T1
GPALN_014649-T1
GPALN_014665-T1
GPALN_014666-T1
GPALN_014669-T1
GPALN_014672-T1
GPALN_014697-T1
GPALN_014697-T2
GPALN_014707-T1
GPALN_014713-T1
GPALN_014715-T1
# name         
GPALN_004921-T1
GPALN_004924-T1
GPALN_004926-T1
GPALN_004930-T1
GPALN_004939-T1
GPALN_004950-T1
GPALN_004951-T1
GPALN_004968-T1
GPALN_004972-T1
GPALN_004974-T1
GPALN_004992-T1
GPALN_005003-T1
GPALN_005008-T1
GPALN_005014-T1
GPALN_005015-T1
GPALN_005016-T1
GPALN_005017-T1
GPALN_005018-T1
GPALN_005021-T1
GPALN_005024-T1
GPALN_005032-T1
GPALN_005038-T1
GPALN_005042-T1
GPALN_005057-T1
GPALN_005058-T1
GPALN_005061-T1
GPALN_005063-T1
GPALN_005067-T1
GPALN_005068-T1
GPALN_005069-T1
GPALN_005070-T1
GPALN_005071-T1
GPALN_005072-T1
GPALN_005074-T1
GPALN_005081-T1
GPALN_005082-T1
GPALN_005088-T1
GPALN_005089-T1
GPALN_005090-T1
GPALN_005097-T1
GPALN_005098-T1
GPALN_005100-T1
GPALN_005101-T1
GPALN_005103-T1
GPALN_005105-T1
GPALN_005106-T1
GPALN_005107-T1
GPALN_005109-T1
GPALN_005113-T1
GPALN_005114-T1
GPALN_005117-T1
GPALN_005120-T1
GPALN_005121-T1
GPALN_005123-T1
GPALN_005129-T1
GPALN_005137-T1
GPALN_005145-T1
GPALN_005147-T1
GPALN_005158-T1
GPALN_005159-T1
GPALN_005160-T1
GPALN_005161-T1
GPALN_005167-T1
GPALN_005172-T1
GPALN_005175-T1
GPALN_005177-T1
GPALN_005180-T1
GPALN_005180-T2
GPALN_005181-T1
GPALN_005185-T1
GPALN_005189-T1
GPALN_005196-T1
GPALN_005198-T1
GPALN_005202-T1
GPALN_005209-T1
GPALN_005209-T2
GPALN_005214-T1
GPALN_005230-T1
GPALN_005234-T1
GPALN_005236-T1
GPALN_005237-T1
GPALN_005238-T1
GPALN_005238-T2
GPALN_005239-T1
GPALN_005242-T1
GPALN_005244-T1
GPALN_005254-T1
GPALN_005263-T1
GPALN_005270-T1
GPALN_005273-T1
GPALN_005276-T1
GPALN_005299-T1
# name         
GPALN_001214-T1
GPALN_001215-T1
GPALN_001216-T1
GPALN_001224-T1
GPALN_001230-T1
GPALN_001251-T1
GPALN_001252-T1
GPALN_001258-T1
GPALN_001281-T1
GPALN_001284-T1
GPALN_001286-T1
GPALN_001298-T1
GPALN_001315-T1
GPALN_001317-T1
GPALN_001318-T1
GPALN_001341-T1
GPALN_001352-T1
GPALN_001355-T1
GPALN_001368-T1
GPALN_001369-T1
GPALN_001449-T1
GPALN_001450-T1
GPALN_001475-T1
GPALN_001478-T1
GPALN_001482-T1
GPALN_001494-T1
GPALN_001496-T1
GPALN_001501-T1
GPALN_001504-T1
GPALN_001505-T1
GPALN_001525-T1
GPALN_001530-T1
GPALN_001534-T1
GPALN_001535-T1
GPALN_001564-T1
GPALN_001590-T1
GPALN_001591-T1
GPALN_001617-T1
GPALN_001625-T1
# name         
GPALN_010593-T1
GPALN_010598-T1
GPALN_010599-T1
GPALN_010600-T1
GPALN_010602-T1
GPALN_010603-T1
GPALN_010613-T1
GPALN_010614-T1
GPALN_010617-T1
GPALN_010620-T1
GPALN_010621-T1
GPALN_010621-T2
GPALN_010625-T1
GPALN_010627-T1
GPALN_010629-T1
GPALN_010632-T1
GPALN_010638-T1
GPALN_010639-T1
GPALN_010643-T1
GPALN_010645-T1
GPALN_010648-T1
GPALN_010650-T1
GPALN_010659-T1
GPALN_010663-T1
GPALN_010664-T1
GPALN_010671-T1
GPALN_010685-T1
GPALN_010702-T1
GPALN_010710-T1
GPALN_010720-T1
GPALN_010723-T1
GPALN_010724-T1
GPALN_010726-T1
GPALN_010727-T1
GPALN_010732-T1
GPALN_010733-T1
GPALN_010737-T1
GPALN_010778-T1
GPALN_010787-T1
GPALN_010789-T1
GPALN_010793-T1
GPALN_010795-T1
GPALN_010798-T1
GPALN_010799-T1
GPALN_010800-T1
GPALN_010801-T1
GPALN_010805-T1
GPALN_010807-T1
GPALN_010809-T1
GPALN_010810-T1
GPALN_010823-T1
GPALN_010824-T1
GPALN_010826-T1
GPALN_010828-T1
GPALN_010829-T1
GPALN_010833-T1
GPALN_010836-T1
GPALN_010861-T1
GPALN_010880-T1
GPALN_010883-T1
GPALN_010884-T1
GPALN_010885-T1
GPALN_010889-T1
GPALN_010895-T1
GPALN_010895-T2
GPALN_010898-T1
GPALN_010903-T1
GPALN_010904-T1
GPALN_010907-T1
GPALN_010908-T1
GPALN_010910-T1
GPALN_010912-T1
GPALN_010945-T1
GPALN_010962-T1
GPALN_010968-T1
GPALN_010970-T1
GPALN_010973-T1
GPALN_010990-T1
GPALN_010991-T1
# name         
GPALN_011426-T1
GPALN_011437-T1
GPALN_011438-T1
GPALN_011444-T1
GPALN_011464-T1
GPALN_011464-T2
GPALN_011495-T1
GPALN_011496-T1
GPALN_011500-T1
GPALN_011506-T1
GPALN_011510-T1
GPALN_011519-T1
GPALN_011522-T1
GPALN_011523-T1
GPALN_011532-T1
GPALN_011535-T1
GPALN_011541-T1
GPALN_011545-T1
GPALN_011546-T1
GPALN_011548-T1
GPALN_011556-T1
GPALN_011558-T1
GPALN_011561-T1
GPALN_011605-T1
GPALN_011607-T1
GPALN_011609-T1
GPALN_011609-T2
GPALN_011610-T1
GPALN_011611-T1
GPALN_011613-T1
GPALN_011614-T1
GPALN_011615-T1
GPALN_011615-T2
GPALN_011615-T3
GPALN_011617-T1
GPALN_011618-T1
GPALN_011620-T1
GPALN_011621-T1
GPALN_011626-T1
GPALN_011627-T1
GPALN_011627-T2
GPALN_011629-T1
GPALN_011634-T1
GPALN_011636-T1
GPALN_011638-T1
GPALN_011649-T1
GPALN_011650-T1
GPALN_011652-T1
GPALN_011693-T1
GPALN_011710-T1
GPALN_011714-T1
GPALN_011715-T1
GPALN_011720-T1
GPALN_011722-T1
GPALN_011724-T1
GPALN_011726-T1
GPALN_011737-T1
GPALN_011738-T1
GPALN_011754-T1
GPALN_011766-T1
GPALN_011775-T1
GPALN_011777-T1
GPALN_011778-T1
GPALN_011797-T1
GPALN_011804-T1
GPALN_011808-T1
GPALN_011812-T1
GPALN_011823-T1
# name         
GPALN_002036-T1
GPALN_002043-T1
GPALN_002052-T1
GPALN_002071-T1
GPALN_002080-T1
GPALN_002082-T1
GPALN_002084-T1
GPALN_002086-T1
GPALN_002098-T1
GPALN_002106-T1
GPALN_002108-T1
GPALN_002117-T1
GPALN_002118-T1
GPALN_002119-T1
GPALN_002125-T1
GPALN_002126-T1
GPALN_002129-T1
GPALN_002129-T2
GPALN_002130-T1
GPALN_002131-T1
GPALN_002133-T1
GPALN_002143-T1
GPALN_002147-T1
GPALN_002148-T1
GPALN_002150-T1
GPALN_002152-T1
GPALN_002162-T1
GPALN_002163-T1
GPALN_002165-T1
GPALN_002168-T1
GPALN_002174-T1
GPALN_002175-T1
GPALN_002176-T1
GPALN_002179-T1
GPALN_002180-T1
GPALN_002181-T1
GPALN_002182-T1
GPALN_002192-T1
GPALN_002194-T1
GPALN_002196-T1
GPALN_002201-T1
GPALN_002204-T1
GPALN_002213-T1
GPALN_002218-T1
GPALN_002219-T1
GPALN_002221-T1
GPALN_002222-T1
GPALN_002223-T1
GPALN_002224-T1
GPALN_002233-T1
GPALN_002246-T1
GPALN_002252-T1
GPALN_002252-T2
GPALN_002252-T3
GPALN_002257-T1
GPALN_002264-T1
GPALN_002267-T1
GPALN_002288-T1
GPALN_002290-T1
GPALN_002294-T1
GPALN_002295-T1
GPALN_002299-T1
GPALN_002300-T1
GPALN_002312-T1
GPALN_002313-T1
GPALN_002314-T1
GPALN_002315-T1
GPALN_002316-T1
GPALN_002318-T1
GPALN_002329-T1
GPALN_002330-T1
GPALN_002334-T1
GPALN_002345-T1
GPALN_002346-T1
GPALN_002348-T1
GPALN_002349-T1
GPALN_002350-T1
GPALN_002352-T1
GPALN_002354-T1
GPALN_002361-T1
GPALN_002364-T1
GPALN_002365-T1
GPALN_002366-T1
GPALN_002370-T1
GPALN_002377-T1
GPALN_002379-T1
GPALN_002382-T1
GPALN_002383-T1
GPALN_002384-T1
GPALN_002387-T1
GPALN_002398-T1
GPALN_002399-T1
GPALN_002418-T1
GPALN_002419-T1
GPALN_002420-T1
GPALN_002422-T1
GPALN_002425-T1
GPALN_002427-T1
GPALN_002432-T1
GPALN_002433-T1
GPALN_002433-T2
# name         
GPALN_011005-T1
GPALN_011018-T1
GPALN_011023-T1
GPALN_011030-T1
GPALN_011036-T1
GPALN_011054-T1
GPALN_011065-T1
GPALN_011069-T1
GPALN_011069-T2
GPALN_011097-T1
GPALN_011124-T1
GPALN_011143-T1
GPALN_011156-T1
GPALN_011159-T1
GPALN_011160-T1
GPALN_011170-T1
GPALN_011172-T1
GPALN_011179-T1
GPALN_011188-T1
GPALN_011197-T1
GPALN_011209-T1
GPALN_011240-T1
GPALN_011260-T1
GPALN_011304-T1
GPALN_011305-T1
GPALN_011309-T1
GPALN_011335-T1
GPALN_011364-T1
GPALN_011375-T1
GPALN_011376-T1
GPALN_011377-T1
GPALN_011378-T1
GPALN_011385-T1
GPALN_011394-T1
GPALN_011396-T1
GPALN_011399-T1
GPALN_011411-T1
GPALN_011411-T2
GPALN_011411-T3
# name         
GPALN_010180-T1
GPALN_010185-T1
GPALN_010199-T1
GPALN_010210-T1
GPALN_010231-T1
GPALN_010232-T1
GPALN_010255-T1
GPALN_010257-T1
GPALN_010270-T1
GPALN_010285-T1
GPALN_010286-T1
GPALN_010289-T1
GPALN_010290-T1
GPALN_010295-T1
GPALN_010296-T1
GPALN_010297-T1
GPALN_010298-T1
GPALN_010299-T1
GPALN_010300-T1
GPALN_010307-T1
GPALN_010316-T1
GPALN_010319-T1
GPALN_010321-T1
GPALN_010330-T1
GPALN_010346-T1
GPALN_010347-T1
GPALN_010348-T1
GPALN_010350-T1
GPALN_010362-T1
GPALN_010388-T1
GPALN_010391-T1
GPALN_010393-T1
GPALN_010403-T1
GPALN_010408-T1
GPALN_010411-T1
GPALN_010413-T1
GPALN_010414-T1
GPALN_010416-T1
GPALN_010417-T1
GPALN_010419-T1
GPALN_010420-T1
GPALN_010424-T1
GPALN_010427-T1
GPALN_010431-T1
GPALN_010432-T1
GPALN_010433-T1
GPALN_010437-T1
GPALN_010444-T1
GPALN_010449-T1
GPALN_010457-T1
GPALN_010458-T1
GPALN_010463-T1
GPALN_010509-T1
GPALN_010510-T1
GPALN_010511-T1
GPALN_010513-T1
GPALN_010517-T1
GPALN_010519-T1
GPALN_010520-T1
GPALN_010524-T1
GPALN_010526-T1
GPALN_010534-T1
GPALN_010535-T1
GPALN_010536-T1
GPALN_010538-T1
GPALN_010540-T1
GPALN_010542-T1
GPALN_010549-T1
GPALN_010552-T1
GPALN_010554-T1
GPALN_010561-T2
GPALN_010568-T1
GPALN_010569-T1
GPALN_010570-T1
GPALN_010571-T1
GPALN_010575-T1
GPALN_010580-T1
GPALN_010582-T1
GPALN_010586-T1
GPALN_010588-T1
# name         
GPALN_016000-T1
GPALN_016014-T1
GPALN_016023-T1
GPALN_016042-T1
GPALN_016050-T1
GPALN_016055-T1
GPALN_016087-T1
GPALN_016087-T2
GPALN_016089-T1
GPALN_016090-T1
GPALN_016091-T1
GPALN_016098-T1
GPALN_016098-T2
GPALN_016098-T3
GPALN_016107-T1
GPALN_016115-T1
GPALN_016116-T1
GPALN_016117-T1
GPALN_016128-T1
GPALN_016131-T1
GPALN_016152-T1
GPALN_016153-T1
GPALN_016157-T1
GPALN_016161-T1
GPALN_016164-T1
GPALN_016166-T1
GPALN_016167-T1
GPALN_016168-T1
GPALN_016172-T1
GPALN_016178-T1
GPALN_016181-T1
GPALN_016182-T1
GPALN_016188-T1
GPALN_016189-T1
GPALN_016191-T1
GPALN_016192-T1
GPALN_016199-T1
GPALN_016201-T1
GPALN_016209-T1
GPALN_016212-T1
GPALN_016213-T1
GPALN_016213-T2
GPALN_016265-T1
GPALN_016268-T1
GPALN_016269-T1
GPALN_016270-T1
GPALN_016277-T1
GPALN_016284-T1
GPALN_016293-T1
GPALN_016297-T1
GPALN_016298-T1
GPALN_016299-T1
GPALN_016330-T1
GPALN_016341-T1
GPALN_016343-T1
GPALN_016347-T1
GPALN_016349-T1
GPALN_016351-T1
GPALN_016359-T1
GPALN_016360-T1
GPALN_016360-T2
GPALN_016373-T1
GPALN_016378-T1
GPALN_016387-T1
# name         
GPALN_015598-T1
GPALN_015600-T1
GPALN_015601-T1
GPALN_015605-T1
GPALN_015616-T1
GPALN_015622-T1
GPALN_015622-T2
GPALN_015629-T1
GPALN_015632-T1
GPALN_015636-T1
GPALN_015640-T1
GPALN_015654-T1
GPALN_015654-T2
GPALN_015655-T1
GPALN_015659-T1
GPALN_015672-T1
GPALN_015675-T1
GPALN_015676-T1
GPALN_015681-T1
GPALN_015686-T1
GPALN_015687-T1
GPALN_015694-T1
GPALN_015696-T1
GPALN_015697-T1
GPALN_015701-T1
GPALN_015704-T1
GPALN_015708-T1
GPALN_015711-T1
GPALN_015712-T1
GPALN_015736-T1
GPALN_015742-T1
GPALN_015743-T1
GPALN_015745-T1
GPALN_015769-T1
GPALN_015771-T1
GPALN_015780-T1
GPALN_015782-T1
GPALN_015799-T1
GPALN_015804-T1
GPALN_015815-T1
GPALN_015819-T1
GPALN_015822-T1
GPALN_015823-T1
GPALN_015826-T1
GPALN_015831-T1
GPALN_015832-T1
GPALN_015833-T1
GPALN_015875-T1
GPALN_015891-T1
GPALN_015910-T1
GPALN_015921-T1
GPALN_015926-T1
GPALN_015927-T1
GPALN_015938-T1
GPALN_015940-T1
GPALN_015947-T1
# name         
GPALN_013060-T1
GPALN_013081-T1
GPALN_013082-T1
GPALN_013083-T1
GPALN_013106-T1
GPALN_013107-T1
GPALN_013108-T1
GPALN_013109-T1
GPALN_013114-T1
GPALN_013115-T1
GPALN_013118-T1
GPALN_013128-T1
GPALN_013133-T1
GPALN_013136-T1
GPALN_013141-T1
GPALN_013142-T1
GPALN_013149-T1
GPALN_013150-T1
GPALN_013158-T1
GPALN_013163-T1
GPALN_013168-T1
GPALN_013170-T1
GPALN_013204-T1
GPALN_013210-T1
GPALN_013214-T1
GPALN_013217-T1
GPALN_013222-T1
GPALN_013224-T1
GPALN_013228-T1
GPALN_013230-T1
GPALN_013233-T1
GPALN_013235-T1
GPALN_013239-T1
GPALN_013242-T1
GPALN_013246-T1
GPALN_013247-T1
GPALN_013248-T1
GPALN_013252-T1
GPALN_013254-T1
GPALN_013256-T1
GPALN_013257-T1
GPALN_013261-T1
GPALN_013263-T1
GPALN_013272-T1
GPALN_013276-T1
GPALN_013277-T1
GPALN_013280-T1
GPALN_013282-T1
GPALN_013283-T1
GPALN_013295-T1
GPALN_013296-T1
GPALN_013297-T1
GPALN_013301-T1
GPALN_013326-T1
GPALN_013347-T1
GPALN_013348-T1
GPALN_013349-T1
GPALN_013350-T1
GPALN_013356-T1
GPALN_013383-T1
GPALN_013384-T1
GPALN_013385-T1
GPALN_013387-T1
GPALN_013399-T1
GPALN_013400-T1
GPALN_013403-T1
GPALN_013403-T2
GPALN_013404-T1
GPALN_013421-T1
GPALN_013423-T1
GPALN_013438-T1
GPALN_013442-T1
GPALN_013448-T1
GPALN_013448-T2
GPALN_013459-T1
GPALN_013463-T1
GPALN_013464-T1
# name         
GPALN_006944-T1
GPALN_006945-T1
GPALN_006949-T1
GPALN_006951-T1
GPALN_006952-T1
GPALN_006963-T1
GPALN_006977-T1
GPALN_006982-T1
GPALN_006988-T1
GPALN_006992-T1
GPALN_006998-T1
GPALN_007023-T1
GPALN_007028-T1
GPALN_007028-T2
GPALN_007033-T1
GPALN_007035-T1
GPALN_007038-T1
GPALN_007039-T1
GPALN_007048-T1
GPALN_007049-T1
GPALN_007051-T1
GPALN_007058-T1
GPALN_007059-T1
GPALN_007060-T1
GPALN_007061-T1
GPALN_007066-T1
GPALN_007070-T1
GPALN_007072-T1
GPALN_007073-T1
GPALN_007074-T1
GPALN_007077-T1
GPALN_007079-T1
GPALN_007082-T1
GPALN_007083-T1
GPALN_007114-T1
GPALN_007115-T1
GPALN_007119-T1
GPALN_007123-T1
GPALN_007129-T1
GPALN_007130-T1
GPALN_007132-T1
GPALN_007142-T1
GPALN_007163-T1
GPALN_007168-T1
GPALN_007170-T1
GPALN_007172-T1
GPALN_007173-T1
GPALN_007176-T1
GPALN_007177-T1
GPALN_007178-T1
GPALN_007179-T1
GPALN_007180-T1
GPALN_007181-T1
GPALN_007185-T1
GPALN_007186-T1
GPALN_007188-T1
GPALN_007189-T1
GPALN_007191-T1
GPALN_007193-T1
GPALN_007194-T1
GPALN_007195-T1
GPALN_007196-T1
GPALN_007199-T1
GPALN_007200-T1
GPALN_007201-T1
GPALN_007202-T1
GPALN_007204-T1
GPALN_007207-T1
GPALN_007207-T2
GPALN_007211-T1
GPALN_007214-T1
GPALN_007215-T1
GPALN_007216-T1
GPALN_007217-T1
GPALN_007218-T1
GPALN_007219-T1
GPALN_007220-T1
GPALN_007221-T1
GPALN_007222-T1
GPALN_007233-T1
GPALN_007237-T1
GPALN_007238-T1
GPALN_007250-T1
GPALN_007259-T1
GPALN_007267-T1
GPALN_007269-T1
GPALN_007286-T1
GPALN_007286-T2
GPALN_007287-T1
GPALN_007293-T1
GPALN_007314-T1
GPALN_007315-T1
GPALN_007317-T1
# name         
GPALN_012254-T1
GPALN_012260-T1
GPALN_012283-T1
GPALN_012284-T1
GPALN_012285-T1
GPALN_012285-T2
GPALN_012287-T1
GPALN_012290-T1
GPALN_012291-T1
GPALN_012293-T1
GPALN_012294-T1
GPALN_012295-T1
GPALN_012296-T1
GPALN_012297-T1
GPALN_012298-T1
GPALN_012299-T1
GPALN_012300-T1
GPALN_012301-T1
GPALN_012302-T1
GPALN_012331-T1
GPALN_012343-T1
GPALN_012357-T1
GPALN_012358-T1
GPALN_012363-T1
GPALN_012366-T1
GPALN_012369-T1
GPALN_012406-T1
GPALN_012409-T1
GPALN_012410-T1
GPALN_012414-T1
GPALN_012415-T1
GPALN_012416-T1
GPALN_012417-T1
GPALN_012435-T1
GPALN_012445-T1
GPALN_012445-T2
GPALN_012452-T1
GPALN_012455-T2
GPALN_012457-T1
GPALN_012465-T1
GPALN_012478-T1
GPALN_012513-T1
GPALN_012523-T1
GPALN_012532-T1
GPALN_012541-T1
GPALN_012550-T1
GPALN_012556-T1
GPALN_012577-T1
GPALN_012580-T1
GPALN_012589-T1
GPALN_012595-T1
GPALN_012602-T1
GPALN_012603-T1
GPALN_012618-T1
GPALN_012619-T1
GPALN_012629-T1
GPALN_012637-T1
GPALN_012641-T1
GPALN_012647-T1
# name         
GPALN_009791-T1
GPALN_009796-T1
GPALN_009815-T1
GPALN_009821-T1
GPALN_009822-T1
GPALN_009823-T1
GPALN_009825-T1
GPALN_009829-T1
GPALN_009835-T1
GPALN_009836-T1
GPALN_009837-T1
GPALN_009839-T1
GPALN_009848-T1
GPALN_009850-T1
GPALN_009860-T1
GPALN_009869-T1
GPALN_009877-T1
GPALN_009879-T1
GPALN_009886-T1
GPALN_009892-T1
GPALN_009893-T1
GPALN_009896-T1
GPALN_009898-T1
GPALN_009899-T1
GPALN_009900-T1
GPALN_009901-T1
GPALN_009902-T1
GPALN_009903-T1
GPALN_009904-T1
GPALN_009905-T1
GPALN_009906-T1
GPALN_009907-T1
GPALN_009908-T1
GPALN_009909-T1
GPALN_009910-T1
GPALN_009911-T1
GPALN_009912-T1
GPALN_009918-T1
GPALN_009944-T1
GPALN_009953-T1
GPALN_009978-T1
GPALN_010012-T1
GPALN_010013-T1
GPALN_010057-T1
GPALN_010059-T1
GPALN_010067-T1
GPALN_010068-T1
GPALN_010073-T1
GPALN_010083-T1
GPALN_010088-T1
GPALN_010093-T1
GPALN_010096-T1
GPALN_010109-T1
GPALN_010112-T1
GPALN_010118-T1
GPALN_010119-T1
GPALN_010120-T1
GPALN_010121-T1
GPALN_010126-T1
GPALN_010127-T1
GPALN_010128-T1
GPALN_010131-T1
GPALN_010136-T1
GPALN_010137-T1
GPALN_010140-T1
GPALN_010148-T1
GPALN_010154-T1
GPALN_010159-T1
GPALN_010160-T1
GPALN_010166-T1
GPALN_010168-T1
GPALN_010171-T1
# name         
GPALN_002860-T1
GPALN_002862-T1
GPALN_002864-T1
GPALN_002879-T1
GPALN_002880-T1
GPALN_002880-T2
GPALN_002881-T1
GPALN_002883-T1
GPALN_002884-T1
GPALN_002885-T1
GPALN_002893-T1
GPALN_002902-T1
GPALN_002905-T1
GPALN_002908-T1
GPALN_002917-T1
GPALN_002941-T1
GPALN_002942-T1
GPALN_002947-T1
GPALN_002955-T1
GPALN_002960-T1
GPALN_002964-T1
GPALN_002965-T1
GPALN_002969-T1
GPALN_002987-T1
GPALN_002988-T1
GPALN_002989-T1
GPALN_002990-T1
GPALN_003006-T1
GPALN_003010-T1
GPALN_003011-T1
GPALN_003015-T1
GPALN_003040-T1
GPALN_003047-T1
GPALN_003055-T1
GPALN_003057-T1
GPALN_003060-T1
GPALN_003061-T1
GPALN_003071-T1
GPALN_003073-T1
GPALN_003075-T1
GPALN_003077-T1
GPALN_003081-T1
GPALN_003083-T1
GPALN_003090-T1
GPALN_003091-T1
GPALN_003092-T1
GPALN_003094-T1
GPALN_003118-T1
GPALN_003175-T1
GPALN_003176-T1
GPALN_003177-T1
GPALN_003178-T1
GPALN_003179-T1
GPALN_003182-T1
GPALN_003185-T1
GPALN_003192-T1
GPALN_003217-T1
GPALN_003218-T1
GPALN_003222-T1
GPALN_003223-T1
GPALN_003246-T1
# name         
GPALN_006522-T1
GPALN_006532-T1
GPALN_006534-T1
GPALN_006555-T1
GPALN_006575-T1
GPALN_006578-T1
GPALN_006581-T1
GPALN_006586-T1
GPALN_006587-T1
GPALN_006588-T1
GPALN_006589-T1
GPALN_006590-T1
GPALN_006594-T1
GPALN_006596-T1
GPALN_006598-T1
GPALN_006599-T1
GPALN_006602-T1
GPALN_006603-T1
GPALN_006604-T1
GPALN_006613-T1
GPALN_006627-T1
GPALN_006629-T1
GPALN_006631-T1
GPALN_006646-T1
GPALN_006650-T1
GPALN_006652-T1
GPALN_006654-T1
GPALN_006656-T1
GPALN_006660-T1
GPALN_006663-T1
GPALN_006687-T1
GPALN_006704-T1
GPALN_006716-T1
GPALN_006719-T1
GPALN_006720-T1
GPALN_006725-T1
GPALN_006726-T1
GPALN_006727-T1
GPALN_006728-T1
GPALN_006730-T1
GPALN_006731-T1
GPALN_006752-T1
GPALN_006754-T1
GPALN_006755-T1
GPALN_006756-T1
GPALN_006759-T1
GPALN_006761-T1
GPALN_006766-T1
GPALN_006769-T1
GPALN_006770-T1
GPALN_006772-T1
GPALN_006775-T1
GPALN_006778-T1
GPALN_006780-T1
GPALN_006782-T1
GPALN_006795-T1
GPALN_006801-T1
GPALN_006802-T1
GPALN_006804-T1
GPALN_006805-T1
GPALN_006806-T1
GPALN_006812-T1
GPALN_006818-T1
GPALN_006819-T1
GPALN_006820-T1
GPALN_006828-T1
GPALN_006839-T1
GPALN_006843-T1
GPALN_006853-T1
GPALN_006854-T1
GPALN_006856-T1
GPALN_006860-T1
GPALN_006864-T1
GPALN_006884-T1
GPALN_006886-T1
GPALN_006907-T1
GPALN_006911-T1
GPALN_006913-T1
GPALN_006914-T1
GPALN_006919-T1
GPALN_006925-T1
GPALN_006932-T1
# name         
GPALN_004487-T1
GPALN_004493-T1
GPALN_004497-T1
GPALN_004506-T1
GPALN_004515-T1
GPALN_004534-T1
GPALN_004552-T1
GPALN_004553-T1
GPALN_004554-T1
GPALN_004555-T1
GPALN_004556-T1
GPALN_004557-T1
GPALN_004559-T1
GPALN_004560-T1
GPALN_004561-T1
GPALN_004564-T1
GPALN_004568-T1
GPALN_004569-T1
GPALN_004571-T1
GPALN_004572-T1
GPALN_004582-T1
GPALN_004585-T1
GPALN_004587-T1
GPALN_004593-T1
GPALN_004595-T1
GPALN_004599-T1
GPALN_004601-T1
GPALN_004602-T1
GPALN_004610-T1
GPALN_004616-T1
GPALN_004626-T1
GPALN_004645-T1
GPALN_004646-T1
GPALN_004647-T1
GPALN_004651-T1
GPALN_004655-T1
GPALN_004665-T1
GPALN_004667-T1
GPALN_004667-T2
GPALN_004668-T1
GPALN_004671-T1
GPALN_004675-T1
GPALN_004676-T1
GPALN_004678-T1
GPALN_004679-T1
GPALN_004681-T1
GPALN_004683-T1
GPALN_004687-T1
GPALN_004689-T1
GPALN_004691-T1
GPALN_004699-T1
GPALN_004700-T1
GPALN_004707-T1
GPALN_004712-T1
GPALN_004729-T1
GPALN_004734-T1
GPALN_004736-T1
GPALN_004738-T1
GPALN_004739-T1
GPALN_004743-T1
GPALN_004749-T1
GPALN_004757-T1
GPALN_004758-T1
GPALN_004770-T1
GPALN_004782-T1
GPALN_004790-T1
GPALN_004794-T1
GPALN_004797-T1
GPALN_004799-T1
GPALN_004800-T1
GPALN_004802-T1
GPALN_004807-T1
GPALN_004810-T1
GPALN_004811-T1
GPALN_004816-T1
GPALN_004818-T1
GPALN_004824-T1
GPALN_004830-T1
GPALN_004840-T1
GPALN_004840-T2
GPALN_004840-T3
GPALN_004842-T1
GPALN_004843-T1
GPALN_004854-T1
GPALN_004856-T1
GPALN_004857-T1
GPALN_004858-T1
GPALN_004862-T1
GPALN_004868-T1
GPALN_004869-T1
GPALN_004876-T1
GPALN_004878-T1
GPALN_004879-T1
GPALN_004880-T1
GPALN_004881-T1
GPALN_004887-T1
GPALN_004897-T1
GPALN_004901-T1
GPALN_004902-T1""".split("\n")

spry = """GPALN_002855-T1
GPALN_002855-T1
GPALN_002855-T1
GPALN_010231-T1
GPALN_010231-T1
GPALN_010231-T1
GPALN_007004-T1
GPALN_007004-T1
GPALN_007004-T1
GPALN_007769-T1
GPALN_007769-T1
GPALN_007771-T1
GPALN_007771-T1
GPALN_010093-T1
GPALN_010093-T1
GPALN_006035-T1
GPALN_006035-T1
GPALN_006035-T1
GPALN_010227-T1
GPALN_010227-T1
GPALN_015852-T1
GPALN_015852-T1
GPALN_007767-T1
GPALN_007767-T1
GPALN_006824-T1
GPALN_006824-T1
GPALN_006830-T1
GPALN_006830-T1
GPALN_006853-T1
GPALN_006853-T1
GPALN_015866-T1
GPALN_015866-T1
GPALN_007766-T1
GPALN_007766-T1
GPALN_015859-T1
GPALN_015859-T1
GPALN_007309-T1
GPALN_007309-T1
GPALN_015115-T1
GPALN_015115-T1
GPALN_015857-T1
GPALN_015857-T1
GPALN_007000-T1
GPALN_007000-T1
GPALN_007013-T1
GPALN_007013-T1
GPALN_007140-T1
GPALN_007140-T1
GPALN_001382-T1
GPALN_001382-T1
GPALN_006597-T1
GPALN_006597-T1
GPALN_009720-T1
GPALN_009720-T1
GPALN_015313-T2
GPALN_015313-T2
GPALN_015313-T1
GPALN_015313-T1
GPALN_010718-T1
GPALN_010718-T1
GPALN_014250-T1
GPALN_014250-T1
GPALN_001749-T1
GPALN_001749-T1
GPALN_012019-T1
GPALN_012019-T1
GPALN_015399-T1
GPALN_015399-T1
GPALN_007015-T1
GPALN_007015-T1
GPALN_011976-T1
GPALN_011976-T1
GPALN_009468-T1
GPALN_009468-T1
GPALN_006981-T1
GPALN_006981-T1
GPALN_015132-T1
GPALN_015132-T1
GPALN_010212-T1
GPALN_010212-T1
GPALN_015314-T1
GPALN_012040-T1
GPALN_012040-T1
GPALN_008170-T1
GPALN_012048-T1
GPALN_012048-T1
GPALN_004633-T1
GPALN_004633-T1
GPALN_011977-T2
GPALN_015845-T1
GPALN_011858-T1
GPALN_015295-T1
GPALN_005620-T1
GPALN_007764-T1
GPALN_015373-T1
GPALN_012367-T1
GPALN_007052-T1
GPALN_007052-T1
GPALN_006596-T1
GPALN_006596-T1
GPALN_012112-T1
GPALN_015632-T1
GPALN_004788-T1
GPALN_014840-T1
GPALN_000148-T1
GPALN_004533-T1
GPALN_011977-T1
GPALN_007026-T1
GPALN_013391-T3
GPALN_015133-T1
GPALN_007153-T1
GPALN_009918-T1
GPALN_013391-T2
GPALN_000146-T1
GPALN_004893-T2
GPALN_010361-T1
GPALN_015628-T1
GPALN_013391-T1
GPALN_003302-T1
GPALN_013350-T1
GPALN_004893-T1
GPALN_006990-T1
GPALN_002300-T1
GPALN_007020-T1
GPALN_007011-T1
GPALN_013171-T1
GPALN_015850-T1
GPALN_007786-T1
GPALN_001353-T1
GPALN_015854-T1
GPALN_015838-T1
GPALN_015122-T1
GPALN_004908-T1
GPALN_015311-T1
GPALN_011963-T1
GPALN_010501-T1
GPALN_010501-T1
GPALN_013393-T1
GPALN_010224-T1
GPALN_013352-T1
GPALN_003214-T1
GPALN_013925-T1
GPALN_007050-T1
GPALN_006837-T1
GPALN_012287-T1
GPALN_006969-T1
GPALN_011999-T1
GPALN_009467-T1
GPALN_015298-T1
GPALN_015315-T1
GPALN_007125-T1
GPALN_007441-T1
GPALN_012194-T1
GPALN_015848-T1
GPALN_009669-T1
GPALN_007146-T1
GPALN_007447-T1
GPALN_007120-T1
GPALN_009870-T1
GPALN_016134-T1
GPALN_009458-T1
GPALN_009917-T1
GPALN_006822-T1
GPALN_015805-T1
GPALN_010063-T1
GPALN_009532-T1
GPALN_006822-T2
GPALN_015309-T1
GPALN_001352-T1
GPALN_011968-T1
GPALN_010097-T1
GPALN_016040-T1
GPALN_007149-T1
GPALN_006856-T1
GPALN_012364-T1
GPALN_006979-T1
GPALN_007012-T1
GPALN_008641-T1
GPALN_012007-T1
GPALN_015301-T1
GPALN_015278-T1
GPALN_015302-T1
GPALN_004455-T1
GPALN_007126-T1
GPALN_012035-T1
GPALN_009815-T1
GPALN_014823-T1
GPALN_007435-T1
GPALN_012056-T1
GPALN_014295-T1
GPALN_014296-T1
GPALN_007778-T1
GPALN_006032-T1
GPALN_013383-T1
GPALN_010392-T1
GPALN_004981-T1
GPALN_007768-T1
GPALN_009709-T1
GPALN_009453-T1
GPALN_014266-T1
GPALN_015013-T1
GPALN_015798-T1
GPALN_007773-T1
GPALN_011765-T1
GPALN_008973-T1
GPALN_011878-T1
GPALN_002323-T1
GPALN_011914-T1
GPALN_005953-T1
GPALN_015569-T1
GPALN_007442-T1
GPALN_011950-T1
GPALN_015101-T1
GPALN_008646-T1
GPALN_012058-T1
GPALN_012060-T1
GPALN_010039-T1
GPALN_016151-T1
GPALN_007128-T1
GPALN_009712-T1
GPALN_011929-T1
GPALN_007777-T1
GPALN_007725-T1
GPALN_015401-T1
GPALN_009714-T1
GPALN_006058-T1
GPALN_013348-T1
GPALN_015858-T1
GPALN_010452-T1
GPALN_009887-T1
GPALN_007439-T1
GPALN_006983-T1
GPALN_006931-T1
GPALN_000147-T1
GPALN_009683-T1
GPALN_007053-T1
GPALN_010793-T1
GPALN_011919-T1
GPALN_003057-T1
GPALN_006635-T1
GPALN_007010-T1
GPALN_015411-T1
GPALN_006991-T1
GPALN_009489-T1
GPALN_015666-T1
GPALN_007131-T1
GPALN_006984-T1
GPALN_009564-T1
GPALN_000144-T1
GPALN_014355-T1
GPALN_007168-T1
GPALN_007452-T1
GPALN_011757-T1
GPALN_013114-T1
GPALN_007790-T1
GPALN_013479-T1
GPALN_010611-T1
GPALN_007132-T1
GPALN_001346-T1
GPALN_004734-T1
GPALN_007129-T1
GPALN_006929-T1
GPALN_015837-T1
GPALN_010789-T1
GPALN_015102-T1
GPALN_007711-T1
GPALN_013074-T1
GPALN_015104-T1
GPALN_007141-T1
GPALN_010597-T1
GPALN_007007-T1
GPALN_010569-T1
GPALN_014357-T1
GPALN_003102-T1
GPALN_012277-T1
GPALN_012020-T1
GPALN_003076-T1
GPALN_007711-T2
GPALN_008217-T1
GPALN_012262-T1
GPALN_007139-T1
GPALN_013379-T1
GPALN_009671-T1
GPALN_007030-T1
GPALN_008654-T1
GPALN_010969-T1
GPALN_007016-T1
GPALN_010791-T1
GPALN_010359-T1
GPALN_002298-T1
GPALN_012111-T1
GPALN_006826-T1
GPALN_010221-T1
GPALN_011759-T1
GPALN_015856-T1
GPALN_015007-T1
GPALN_015274-T1
GPALN_015405-T2
GPALN_007774-T1
GPALN_015398-T1
GPALN_015405-T1
GPALN_003103-T1
GPALN_015868-T1
GPALN_013353-T1
GPALN_010653-T1
GPALN_010222-T1
GPALN_011827-T1
GPALN_010211-T1
GPALN_009672-T1
GPALN_011794-T1
GPALN_014373-T1
GPALN_006861-T1
GPALN_010970-T1
GPALN_014256-T1
GPALN_014477-T1
GPALN_012306-T1
GPALN_014344-T1
GPALN_010441-T1
GPALN_015107-T1
GPALN_006775-T1
GPALN_009476-T1
GPALN_015407-T1
GPALN_013385-T1
GPALN_014397-T1
GPALN_007454-T1
GPALN_010208-T1
GPALN_007491-T1
GPALN_007458-T1
GPALN_014398-T1
GPALN_014257-T1
GPALN_004906-T2
GPALN_013168-T1
GPALN_015830-T1
GPALN_007009-T1
GPALN_006939-T1
GPALN_003112-T1
GPALN_011948-T1
GPALN_009445-T1
GPALN_004891-T1
GPALN_011816-T1
GPALN_007019-T1
GPALN_012392-T1
GPALN_004877-T1
GPALN_007421-T1
GPALN_002297-T1
GPALN_006828-T1
GPALN_004906-T1
GPALN_015377-T1
GPALN_002290-T1
GPALN_015087-T1
GPALN_007780-T1
GPALN_015803-T1
GPALN_004761-T1
GPALN_004918-T1
GPALN_013437-T1
GPALN_014395-T1
GPALN_013127-T1
GPALN_014233-T1
GPALN_009915-T1
GPALN_010232-T1
GPALN_013436-T1
GPALN_013111-T1
GPALN_003794-T1
GPALN_006818-T1
GPALN_002293-T1
GPALN_004265-T1
GPALN_007006-T1
GPALN_010645-T1
GPALN_004783-T1
GPALN_004471-T1
GPALN_004003-T1
GPALN_009871-T1
GPALN_006968-T1
GPALN_011761-T1
GPALN_012064-T1
GPALN_004632-T1
GPALN_007045-T2
GPALN_006836-T1
GPALN_007045-T1
GPALN_015814-T1
GPALN_006839-T1
GPALN_010216-T1
GPALN_004313-T1
GPALN_014246-T1
GPALN_011823-T1
GPALN_015800-T1
GPALN_004310-T1
GPALN_009924-T1
GPALN_014271-T1
GPALN_015863-T1
GPALN_004900-T1
GPALN_007024-T1
GPALN_004627-T1
GPALN_012028-T1
GPALN_007047-T2
GPALN_007144-T1
GPALN_014267-T1
GPALN_012062-T1
GPALN_010085-T1
GPALN_007047-T1
GPALN_007459-T1
GPALN_005857-T1
GPALN_006975-T1
GPALN_013342-T1
GPALN_009446-T1
GPALN_009459-T1
GPALN_006664-T1
GPALN_015173-T1
GPALN_003793-T1
GPALN_016362-T1
GPALN_014864-T1
GPALN_011962-T1
GPALN_007187-T1
GPALN_009465-T1
GPALN_009681-T1
GPALN_001077-T1
GPALN_014560-T1
GPALN_015109-T1
GPALN_006036-T1
GPALN_010693-T1
GPALN_007155-T1
GPALN_011811-T1
GPALN_015810-T1
GPALN_014963-T1
GPALN_003304-T1
GPALN_013480-T1
GPALN_014548-T1
GPALN_012985-T1
GPALN_007122-T1
GPALN_006033-T1
GPALN_014574-T1
GPALN_009811-T1
GPALN_012365-T1
GPALN_007894-T1
GPALN_009806-T1
GPALN_002979-T1
GPALN_010967-T1
GPALN_001750-T1
GPALN_015195-T1
GPALN_009814-T1
GPALN_009697-T1
GPALN_007424-T1
GPALN_014249-T1
GPALN_006858-T1
GPALN_014247-T1
GPALN_009673-T1
GPALN_011815-T1
GPALN_004735-T1
GPALN_010366-T1
GPALN_006831-T1
GPALN_002288-T1
GPALN_016160-T1
GPALN_014762-T1
GPALN_016362-T2
GPALN_004809-T1
GPALN_006840-T1
GPALN_015124-T1
GPALN_007002-T1
GPALN_004920-T1
GPALN_010236-T1
GPALN_007445-T1
GPALN_006940-T1
GPALN_009473-T1
GPALN_007150-T1
GPALN_007025-T1
GPALN_011821-T1
GPALN_011826-T1
GPALN_009462-T1
GPALN_015375-T1
GPALN_012006-T1
GPALN_007056-T1
GPALN_014573-T1
GPALN_013178-T1
GPALN_006851-T1
GPALN_002277-T1
GPALN_014347-T1
GPALN_007136-T1
GPALN_004910-T1
GPALN_015280-T1
GPALN_011937-T1
GPALN_011938-T1
GPALN_006860-T1
GPALN_013343-T1
GPALN_010968-T1
GPALN_004897-T1
GPALN_015813-T1
GPALN_006842-T1
GPALN_006845-T1
GPALN_002281-T1
GPALN_015371-T1
GPALN_003114-T1
GPALN_015836-T1
GPALN_012018-T1
GPALN_004881-T1
GPALN_014570-T1
GPALN_011631-T1
GPALN_015867-T1
GPALN_015189-T1
GPALN_010261-T1
GPALN_015108-T1
GPALN_013414-T1
GPALN_013346-T1
GPALN_015368-T1
GPALN_015362-T1
GPALN_004896-T1
GPALN_001780-T1
GPALN_004635-T1
GPALN_015111-T1
GPALN_011695-T1
GPALN_006833-T1
GPALN_007451-T1
GPALN_009472-T1
GPALN_013339-T1
GPALN_007456-T1
GPALN_009808-T1""".split("\n")

secreted_name = set([])

for i in phobius:
    name = i.split()[0]
    name = name.rstrip()
    tm = i.split()[1]
    if name.endswith("-T1"):
        if tm == "0":
            secreted_name.add(name)

secreted_name_sigP = set([])

for i in signalP:
    if i.startswith("#"):
        continue
    name = i
    name = name.rstrip()
    if name.endswith("-T1"):
        secreted_name_sigP.add(name)

spry_name = set([])
for i in spry:
    name = i.rstrip()
    if name.endswith("-T1"):
        spry_name.add(name)

sig_and_phobius = secreted_name.union(secreted_name_sigP)
SPRYSECS = spry_name.intersection(sig_and_phobius)
print(len(SPRYSECS))

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:


$ python name of scrip -h

make sure you run this in the folder where the file is...


"""

parser = OptionParser(usage=usage)
parser.add_option("-i", dest="infile", default="all_putative_effectors.txt",
                    help="tab output from scanMotifs.pl")
                  
parser.add_option("-f", dest="fasta", 
                    default="./data/Gpal_newton_annotated.AA.fasta",
                    help="fasta with annotation in the desctption")
                  
parser.add_option("-o", "--out", dest="out", 
                default="results_annotations.txt",
                help="Output filename",
                metavar="FILE")
                  
parser.add_option("-l", "--log", dest="log", 
                default="logger.out",
                help="Output filename",
                metavar="FILE")
                  
parser.add_option("-c", "--caz", dest="caz", 
                 default="./data/annotations.dbCAN.txt",
                   help="Output filename",
                metavar="FILE")

parser.add_option("-d", "--dup", dest="duplication", 
                    default="./data/gp_gp.gene_type",
                    help="duplication file",
                    metavar="FILE")

parser.add_option("--dpi14_vs_Gp_dpi21_Gp_dpi14_UP", dest="dpi14_vs_Gp_dpi21_Gp_dpi14_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_dpi21.edgeR.DE_results.P1e-3_C2_Gp_dpi14_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi14_vs_Gp_dpi21_Gp_dpi21_UP", dest="dpi14_vs_Gp_dpi21_Gp_dpi21_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_dpi21.edgeR.DE_results.P1e-3_C2_Gp_dpi21_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi14_vs_Gp_dpi28_Gp_dpi14_UP", dest="dpi14_vs_Gp_dpi28_Gp_dpi14_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_dpi28.edgeR.DE_results.P1e-3_C2_Gp_dpi14_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi14_vs_Gp_dpi28_Gp_dpi28_UP", dest="dpi14_vs_Gp_dpi28_Gp_dpi28_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_dpi28.edgeR.DE_results.P1e-3_C2_Gp_dpi28_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi14_vs_Gp_dpi35_Gp_dpi14_UP", dest="dpi14_vs_Gp_dpi35_Gp_dpi14_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_dpi35.edgeR.DE_results.P1e-3_C2_Gp_dpi14_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi14_vs_Gp_dpi35_Gp_dpi35_UP", dest="dpi14_vs_Gp_dpi35_Gp_dpi35_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_dpi35.edgeR.DE_results.P1e-3_C2_Gp_dpi35_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi14_vs_Gp_dpi7_Gp_dpi14_UP", dest="dpi14_vs_Gp_dpi7_Gp_dpi14_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_dpi7.edgeR.DE_results.P1e-3_C2_Gp_dpi14_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi14_vs_Gp_dpi7_Gp_dpi7_UP", dest="dpi14_vs_Gp_dpi7_Gp_dpi7_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_dpi7.edgeR.DE_results.P1e-3_C2_Gp_dpi7_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi14_vs_Gp_EGG_Gp_dpi14_UP", dest="dpi14_vs_Gp_EGG_Gp_dpi14_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_EGG.edgeR.DE_results.P1e-3_C2_Gp_dpi14_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi14_vs_Gp_EGG_Gp_EGG_UP", dest="dpi14_vs_Gp_EGG_Gp_EGG_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_EGG.edgeR.DE_results.P1e-3_C2_Gp_EGG_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi14_vs_Gp_J2_Gp_dpi14_UP", dest="dpi14_vs_Gp_J2_Gp_dpi14_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_J2.edgeR.DE_results.P1e-3_C2_Gp_dpi14_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi14_vs_Gp_J2_Gp_J2_UP", dest="dpi14_vs_Gp_J2_Gp_J2_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_J2.edgeR.DE_results.P1e-3_C2_Gp_J2_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi14_vs_Gp_MALE_Gp_dpi14_UP", dest="dpi14_vs_Gp_MALE_Gp_dpi14_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_dpi14_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi14_vs_Gp_MALE_Gp_MALE_UP", dest="dpi14_vs_Gp_MALE_Gp_MALE_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi14_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_MALE_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi21_vs_Gp_dpi28_Gp_dpi21_UP", dest="dpi21_vs_Gp_dpi28_Gp_dpi21_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi21_vs_Gp_dpi28.edgeR.DE_results.P1e-3_C2_Gp_dpi21_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi21_vs_Gp_dpi28_Gp_dpi28_UP", dest="dpi21_vs_Gp_dpi28_Gp_dpi28_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi21_vs_Gp_dpi28.edgeR.DE_results.P1e-3_C2_Gp_dpi28_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi21_vs_Gp_dpi35_Gp_dpi21_UP", dest="dpi21_vs_Gp_dpi35_Gp_dpi21_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi21_vs_Gp_dpi35.edgeR.DE_results.P1e-3_C2_Gp_dpi21_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi21_vs_Gp_dpi35_Gp_dpi35_UP", dest="dpi21_vs_Gp_dpi35_Gp_dpi35_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi21_vs_Gp_dpi35.edgeR.DE_results.P1e-3_C2_Gp_dpi35_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi21_vs_Gp_dpi7_Gp_dpi21_UP", dest="dpi21_vs_Gp_dpi7_Gp_dpi21_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi21_vs_Gp_dpi7.edgeR.DE_results.P1e-3_C2_Gp_dpi21_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi21_vs_Gp_dpi7_Gp_dpi7_UP", dest="dpi21_vs_Gp_dpi7_Gp_dpi7_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi21_vs_Gp_dpi7.edgeR.DE_results.P1e-3_C2_Gp_dpi7_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi21_vs_Gp_EGG_Gp_dpi21_UP", dest="dpi21_vs_Gp_EGG_Gp_dpi21_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi21_vs_Gp_EGG.edgeR.DE_results.P1e-3_C2_Gp_dpi21_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi21_vs_Gp_EGG_Gp_EGG_UP", dest="dpi21_vs_Gp_EGG_Gp_EGG_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi21_vs_Gp_EGG.edgeR.DE_results.P1e-3_C2_Gp_EGG_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi21_vs_Gp_J2_Gp_dpi21_UP", dest="dpi21_vs_Gp_J2_Gp_dpi21_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi21_vs_Gp_J2.edgeR.DE_results.P1e-3_C2_Gp_dpi21_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi21_vs_Gp_J2_Gp_J2_UP", dest="dpi21_vs_Gp_J2_Gp_J2_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi21_vs_Gp_J2.edgeR.DE_results.P1e-3_C2_Gp_J2_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi21_vs_Gp_MALE_Gp_dpi21_UP", dest="dpi21_vs_Gp_MALE_Gp_dpi21_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi21_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_dpi21_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi21_vs_Gp_MALE_Gp_MALE_UP", dest="dpi21_vs_Gp_MALE_Gp_MALE_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi21_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_MALE_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi28_vs_Gp_dpi35_Gp_dpi28_UP", dest="dpi28_vs_Gp_dpi35_Gp_dpi28_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi28_vs_Gp_dpi35.edgeR.DE_results.P1e-3_C2_Gp_dpi28_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi28_vs_Gp_dpi35_Gp_dpi35_UP", dest="dpi28_vs_Gp_dpi35_Gp_dpi35_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi28_vs_Gp_dpi35.edgeR.DE_results.P1e-3_C2_Gp_dpi35_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi28_vs_Gp_dpi7_Gp_dpi28_UP", dest="dpi28_vs_Gp_dpi7_Gp_dpi28_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi28_vs_Gp_dpi7.edgeR.DE_results.P1e-3_C2_Gp_dpi28_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi28_vs_Gp_dpi7_Gp_dpi7_UP", dest="dpi28_vs_Gp_dpi7_Gp_dpi7_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi28_vs_Gp_dpi7.edgeR.DE_results.P1e-3_C2_Gp_dpi7_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi28_vs_Gp_EGG_Gp_dpi28_UP", dest="dpi28_vs_Gp_EGG_Gp_dpi28_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi28_vs_Gp_EGG.edgeR.DE_results.P1e-3_C2_Gp_dpi28_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi28_vs_Gp_EGG_Gp_EGG_UP", dest="dpi28_vs_Gp_EGG_Gp_EGG_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi28_vs_Gp_EGG.edgeR.DE_results.P1e-3_C2_Gp_EGG_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi28_vs_Gp_J2_Gp_dpi28_UP", dest="dpi28_vs_Gp_J2_Gp_dpi28_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi28_vs_Gp_J2.edgeR.DE_results.P1e-3_C2_Gp_dpi28_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi28_vs_Gp_J2_Gp_J2_UP", dest="dpi28_vs_Gp_J2_Gp_J2_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi28_vs_Gp_J2.edgeR.DE_results.P1e-3_C2_Gp_J2_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi28_vs_Gp_MALE_Gp_dpi28_UP", dest="dpi28_vs_Gp_MALE_Gp_dpi28_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi28_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_dpi28_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi28_vs_Gp_MALE_Gp_MALE_UP", dest="dpi28_vs_Gp_MALE_Gp_MALE_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi28_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_MALE_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi35_vs_Gp_dpi7_Gp_dpi35_UP", dest="dpi35_vs_Gp_dpi7_Gp_dpi35_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi35_vs_Gp_dpi7.edgeR.DE_results.P1e-3_C2_Gp_dpi35_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi35_vs_Gp_dpi7_Gp_dpi7_UP", dest="dpi35_vs_Gp_dpi7_Gp_dpi7_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi35_vs_Gp_dpi7.edgeR.DE_results.P1e-3_C2_Gp_dpi7_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi35_vs_Gp_EGG_Gp_dpi35_UP", dest="dpi35_vs_Gp_EGG_Gp_dpi35_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi35_vs_Gp_EGG.edgeR.DE_results.P1e-3_C2_Gp_dpi35_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi35_vs_Gp_EGG_Gp_EGG_UP", dest="dpi35_vs_Gp_EGG_Gp_EGG_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi35_vs_Gp_EGG.edgeR.DE_results.P1e-3_C2_Gp_EGG_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi35_vs_Gp_J2_Gp_dpi35_UP", dest="dpi35_vs_Gp_J2_Gp_dpi35_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi35_vs_Gp_J2.edgeR.DE_results.P1e-3_C2_Gp_dpi35_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi35_vs_Gp_J2_Gp_J2_UP", dest="dpi35_vs_Gp_J2_Gp_J2_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi35_vs_Gp_J2.edgeR.DE_results.P1e-3_C2_Gp_J2_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi35_vs_Gp_MALE_Gp_dpi35_UP", dest="dpi35_vs_Gp_MALE_Gp_dpi35_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi35_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_dpi35_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi35_vs_Gp_MALE_Gp_MALE_UP", dest="dpi35_vs_Gp_MALE_Gp_MALE_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi35_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_MALE_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi7_vs_Gp_EGG_Gp_dpi7_UP", dest="dpi7_vs_Gp_EGG_Gp_dpi7_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi7_vs_Gp_EGG.edgeR.DE_results.P1e-3_C2_Gp_dpi7_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi7_vs_Gp_EGG_Gp_EGG_UP", dest="dpi7_vs_Gp_EGG_Gp_EGG_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi7_vs_Gp_EGG.edgeR.DE_results.P1e-3_C2_Gp_EGG_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi7_vs_Gp_J2_Gp_dpi7_UP", dest="dpi7_vs_Gp_J2_Gp_dpi7_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi7_vs_Gp_J2.edgeR.DE_results.P1e-3_C2_Gp_dpi7_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi7_vs_Gp_J2_Gp_J2_UP", dest="dpi7_vs_Gp_J2_Gp_J2_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi7_vs_Gp_J2.edgeR.DE_results.P1e-3_C2_Gp_J2_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi7_vs_Gp_MALE_Gp_dpi7_UP", dest="dpi7_vs_Gp_MALE_Gp_dpi7_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi7_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_dpi7_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--dpi7_vs_Gp_MALE_Gp_MALE_UP", dest="dpi7_vs_Gp_MALE_Gp_MALE_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_dpi7_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_MALE_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--EGG_vs_Gp_J2_Gp_EGG_UP", dest="EGG_vs_Gp_J2_Gp_EGG_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_EGG_vs_Gp_J2.edgeR.DE_results.P1e-3_C2_Gp_EGG_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--EGG_vs_Gp_J2_Gp_J2_UP", dest="EGG_vs_Gp_J2_Gp_J2_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_EGG_vs_Gp_J2.edgeR.DE_results.P1e-3_C2_Gp_J2_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--EGG_vs_Gp_MALE_Gp_EGG_UP", dest="EGG_vs_Gp_MALE_Gp_EGG_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_EGG_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_EGG_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--EGG_vs_Gp_MALE_Gp_MALE_UP", dest="EGG_vs_Gp_MALE_Gp_MALE_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_EGG_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_MALE_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--J2_vs_Gp_MALE_Gp_J2_UP", dest="J2_vs_Gp_MALE_Gp_J2_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_J2_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_J2_UP.subset",  
                  help="DE file",
                  metavar="FILE")
 
parser.add_option("--J2_vs_Gp_MALE_Gp_MALE_UP", dest="J2_vs_Gp_MALE_Gp_MALE_UP", 
                  default="./DE/Gp_genes.counts.matrix_Gp_J2_vs_Gp_MALE.edgeR.DE_results.P1e-3_C2_Gp_MALE_UP.subset",  
                  help="DE file",
                  metavar="FILE")
                  
                  


(options, args) = parser.parse_args()
 
log_out = options.log
 
# Run as script
if __name__ == '__main__':
    # Set up logging
    logger = logging.getLogger('parse_scan_motifs.py: %s' % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    try:
        logging.basicConfig(filename='example.log',level=logging.DEBUG)
        logstream = open(log_out, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        logger.error("Could not open %s for logging", logger)
        sys.exit(1)
    # Report input arguments
    logger.info(sys.version_info)
    logger.info("Command-line: %s", ' '.join(sys.argv))
    logger.info("Starting testing: %s", time.asctime())
    if not os.path.isfile(options.fasta):
        logger.warning("Input genome file not found: %s" % options.fasta)
        sys.exit("Input genome file not found: %s" % options.fasta)
    if not os.path.isfile(options.infile):
        logger.warning("Input BAM file not found: %s" % options.infile)
        sys.exit("Input BAM file not found: %s" % options.infile)
    # gene the gene annota as a dict
    gene_annot = index_annotation(options.fasta)
    # count the number of motifs per gene
    gene_motif_counter = parse_file(options.infile)
    # get th GH domian info
    gene_GH_domain = parse_GH_file(options.caz)
    ###################################################################
    ###################################################################
    dpi14_vs_Gp_dpi21_Gp_dpi14_UP = parse_DE_file(options.dpi14_vs_Gp_dpi21_Gp_dpi14_UP)
    dpi14_vs_Gp_dpi21_Gp_dpi21_UP = parse_DE_file(options.dpi14_vs_Gp_dpi21_Gp_dpi21_UP)
    dpi14_vs_Gp_dpi28_Gp_dpi14_UP = parse_DE_file(options.dpi14_vs_Gp_dpi28_Gp_dpi14_UP)
    dpi14_vs_Gp_dpi28_Gp_dpi28_UP = parse_DE_file(options.dpi14_vs_Gp_dpi28_Gp_dpi28_UP)
    dpi14_vs_Gp_dpi35_Gp_dpi14_UP = parse_DE_file(options.dpi14_vs_Gp_dpi35_Gp_dpi14_UP)
    dpi14_vs_Gp_dpi35_Gp_dpi35_UP = parse_DE_file(options.dpi14_vs_Gp_dpi35_Gp_dpi35_UP)
    dpi14_vs_Gp_dpi7_Gp_dpi14_UP = parse_DE_file(options.dpi14_vs_Gp_dpi7_Gp_dpi14_UP)
    dpi14_vs_Gp_dpi7_Gp_dpi7_UP = parse_DE_file(options.dpi14_vs_Gp_dpi7_Gp_dpi7_UP)
    dpi14_vs_Gp_EGG_Gp_dpi14_UP = parse_DE_file(options.dpi14_vs_Gp_EGG_Gp_dpi14_UP)
    dpi14_vs_Gp_EGG_Gp_EGG_UP = parse_DE_file(options.dpi14_vs_Gp_EGG_Gp_EGG_UP)
    dpi14_vs_Gp_J2_Gp_dpi14_UP = parse_DE_file(options.dpi14_vs_Gp_J2_Gp_dpi14_UP)
    dpi14_vs_Gp_J2_Gp_J2_UP = parse_DE_file(options.dpi14_vs_Gp_J2_Gp_J2_UP)
    dpi14_vs_Gp_MALE_Gp_dpi14_UP = parse_DE_file(options.dpi14_vs_Gp_MALE_Gp_dpi14_UP)
    dpi14_vs_Gp_MALE_Gp_MALE_UP = parse_DE_file(options.dpi14_vs_Gp_MALE_Gp_MALE_UP)
    dpi21_vs_Gp_dpi28_Gp_dpi21_UP = parse_DE_file(options.dpi21_vs_Gp_dpi28_Gp_dpi21_UP)
    dpi21_vs_Gp_dpi28_Gp_dpi28_UP = parse_DE_file(options.dpi21_vs_Gp_dpi28_Gp_dpi28_UP)
    dpi21_vs_Gp_dpi35_Gp_dpi21_UP = parse_DE_file(options.dpi21_vs_Gp_dpi35_Gp_dpi21_UP)
    dpi21_vs_Gp_dpi35_Gp_dpi35_UP = parse_DE_file(options.dpi21_vs_Gp_dpi35_Gp_dpi35_UP)
    dpi21_vs_Gp_dpi7_Gp_dpi21_UP = parse_DE_file(options.dpi21_vs_Gp_dpi7_Gp_dpi21_UP)
    dpi21_vs_Gp_dpi7_Gp_dpi7_UP = parse_DE_file(options.dpi21_vs_Gp_dpi7_Gp_dpi7_UP)
    dpi21_vs_Gp_EGG_Gp_dpi21_UP = parse_DE_file(options.dpi21_vs_Gp_EGG_Gp_dpi21_UP)
    dpi21_vs_Gp_EGG_Gp_EGG_UP = parse_DE_file(options.dpi21_vs_Gp_EGG_Gp_EGG_UP)
    dpi21_vs_Gp_J2_Gp_dpi21_UP = parse_DE_file(options.dpi21_vs_Gp_J2_Gp_dpi21_UP)
    dpi21_vs_Gp_J2_Gp_J2_UP = parse_DE_file(options.dpi21_vs_Gp_J2_Gp_J2_UP)
    dpi21_vs_Gp_MALE_Gp_dpi21_UP = parse_DE_file(options.dpi21_vs_Gp_MALE_Gp_dpi21_UP)
    dpi21_vs_Gp_MALE_Gp_MALE_UP = parse_DE_file(options.dpi21_vs_Gp_MALE_Gp_MALE_UP)
    dpi28_vs_Gp_dpi35_Gp_dpi28_UP = parse_DE_file(options.dpi28_vs_Gp_dpi35_Gp_dpi28_UP)
    dpi28_vs_Gp_dpi35_Gp_dpi35_UP = parse_DE_file(options.dpi28_vs_Gp_dpi35_Gp_dpi35_UP)
    dpi28_vs_Gp_dpi7_Gp_dpi28_UP = parse_DE_file(options.dpi28_vs_Gp_dpi7_Gp_dpi28_UP)
    dpi28_vs_Gp_dpi7_Gp_dpi7_UP = parse_DE_file(options.dpi28_vs_Gp_dpi7_Gp_dpi7_UP)
    dpi28_vs_Gp_EGG_Gp_dpi28_UP = parse_DE_file(options.dpi28_vs_Gp_EGG_Gp_dpi28_UP)
    dpi28_vs_Gp_EGG_Gp_EGG_UP = parse_DE_file(options.dpi28_vs_Gp_EGG_Gp_EGG_UP)
    dpi28_vs_Gp_J2_Gp_dpi28_UP = parse_DE_file(options.dpi28_vs_Gp_J2_Gp_dpi28_UP)
    dpi28_vs_Gp_J2_Gp_J2_UP = parse_DE_file(options.dpi28_vs_Gp_J2_Gp_J2_UP)
    dpi28_vs_Gp_MALE_Gp_dpi28_UP = parse_DE_file(options.dpi28_vs_Gp_MALE_Gp_dpi28_UP)
    dpi28_vs_Gp_MALE_Gp_MALE_UP = parse_DE_file(options.dpi28_vs_Gp_MALE_Gp_MALE_UP)
    dpi35_vs_Gp_dpi7_Gp_dpi35_UP = parse_DE_file(options.dpi35_vs_Gp_dpi7_Gp_dpi35_UP)
    dpi35_vs_Gp_dpi7_Gp_dpi7_UP = parse_DE_file(options.dpi35_vs_Gp_dpi7_Gp_dpi7_UP)
    dpi35_vs_Gp_EGG_Gp_dpi35_UP = parse_DE_file(options.dpi35_vs_Gp_EGG_Gp_dpi35_UP)
    dpi35_vs_Gp_EGG_Gp_EGG_UP = parse_DE_file(options.dpi35_vs_Gp_EGG_Gp_EGG_UP)
    dpi35_vs_Gp_J2_Gp_dpi35_UP = parse_DE_file(options.dpi35_vs_Gp_J2_Gp_dpi35_UP)
    dpi35_vs_Gp_J2_Gp_J2_UP = parse_DE_file(options.dpi35_vs_Gp_J2_Gp_J2_UP)
    dpi35_vs_Gp_MALE_Gp_dpi35_UP = parse_DE_file(options.dpi35_vs_Gp_MALE_Gp_dpi35_UP)
    dpi35_vs_Gp_MALE_Gp_MALE_UP = parse_DE_file(options.dpi35_vs_Gp_MALE_Gp_MALE_UP)
    dpi7_vs_Gp_EGG_Gp_dpi7_UP = parse_DE_file(options.dpi7_vs_Gp_EGG_Gp_dpi7_UP)
    dpi7_vs_Gp_EGG_Gp_EGG_UP = parse_DE_file(options.dpi7_vs_Gp_EGG_Gp_EGG_UP)
    dpi7_vs_Gp_J2_Gp_dpi7_UP = parse_DE_file(options.dpi7_vs_Gp_J2_Gp_dpi7_UP)
    dpi7_vs_Gp_J2_Gp_J2_UP = parse_DE_file(options.dpi7_vs_Gp_J2_Gp_J2_UP)
    dpi7_vs_Gp_MALE_Gp_dpi7_UP = parse_DE_file(options.dpi7_vs_Gp_MALE_Gp_dpi7_UP)
    dpi7_vs_Gp_MALE_Gp_MALE_UP = parse_DE_file(options.dpi7_vs_Gp_MALE_Gp_MALE_UP)
    EGG_vs_Gp_J2_Gp_EGG_UP = parse_DE_file(options.EGG_vs_Gp_J2_Gp_EGG_UP)
    EGG_vs_Gp_J2_Gp_J2_UP = parse_DE_file(options.EGG_vs_Gp_J2_Gp_J2_UP)
    EGG_vs_Gp_MALE_Gp_EGG_UP = parse_DE_file(options.EGG_vs_Gp_MALE_Gp_EGG_UP)
    EGG_vs_Gp_MALE_Gp_MALE_UP = parse_DE_file(options.EGG_vs_Gp_MALE_Gp_MALE_UP)
    J2_vs_Gp_MALE_Gp_J2_UP = parse_DE_file(options.J2_vs_Gp_MALE_Gp_J2_UP)
    J2_vs_Gp_MALE_Gp_MALE_UP = parse_DE_file(options.J2_vs_Gp_MALE_Gp_MALE_UP)

    
    # get the ducpliation types
    gene_duplication = parse_duplication_file(options.duplication)
    # print out the data to a file
    print_out_data(gene_annot, gene_motif_counter, 
                   SPRYSECS, sig_and_phobius, 
                   gene_GH_domain, gene_duplication,
                   dpi14_vs_Gp_dpi21_Gp_dpi14_UP, dpi14_vs_Gp_dpi21_Gp_dpi21_UP, 
                   dpi14_vs_Gp_dpi28_Gp_dpi14_UP, dpi14_vs_Gp_dpi28_Gp_dpi28_UP, 
                   dpi14_vs_Gp_dpi35_Gp_dpi14_UP, dpi14_vs_Gp_dpi35_Gp_dpi35_UP, 
                   dpi14_vs_Gp_dpi7_Gp_dpi14_UP, dpi14_vs_Gp_dpi7_Gp_dpi7_UP, 
                   dpi14_vs_Gp_EGG_Gp_dpi14_UP, dpi14_vs_Gp_EGG_Gp_EGG_UP, 
                   dpi14_vs_Gp_J2_Gp_dpi14_UP, dpi14_vs_Gp_J2_Gp_J2_UP, 
                   dpi14_vs_Gp_MALE_Gp_dpi14_UP, dpi14_vs_Gp_MALE_Gp_MALE_UP, 
                   dpi21_vs_Gp_dpi28_Gp_dpi21_UP, dpi21_vs_Gp_dpi28_Gp_dpi28_UP, 
                   dpi21_vs_Gp_dpi35_Gp_dpi21_UP, dpi21_vs_Gp_dpi35_Gp_dpi35_UP,
                   dpi21_vs_Gp_dpi7_Gp_dpi21_UP, dpi21_vs_Gp_dpi7_Gp_dpi7_UP, 
                   dpi21_vs_Gp_EGG_Gp_dpi21_UP, dpi21_vs_Gp_EGG_Gp_EGG_UP, 
                   dpi21_vs_Gp_J2_Gp_dpi21_UP, dpi21_vs_Gp_J2_Gp_J2_UP, 
                   dpi21_vs_Gp_MALE_Gp_dpi21_UP, dpi21_vs_Gp_MALE_Gp_MALE_UP, 
                   dpi28_vs_Gp_dpi35_Gp_dpi28_UP, dpi28_vs_Gp_dpi35_Gp_dpi35_UP, 
                   dpi28_vs_Gp_dpi7_Gp_dpi28_UP, dpi28_vs_Gp_dpi7_Gp_dpi7_UP, 
                   dpi28_vs_Gp_EGG_Gp_dpi28_UP, dpi28_vs_Gp_EGG_Gp_EGG_UP, 
                   dpi28_vs_Gp_J2_Gp_dpi28_UP, dpi28_vs_Gp_J2_Gp_J2_UP, 
                   dpi28_vs_Gp_MALE_Gp_dpi28_UP, dpi28_vs_Gp_MALE_Gp_MALE_UP, 
                   dpi35_vs_Gp_dpi7_Gp_dpi35_UP, dpi35_vs_Gp_dpi7_Gp_dpi7_UP, 
                   dpi35_vs_Gp_EGG_Gp_dpi35_UP, dpi35_vs_Gp_EGG_Gp_EGG_UP, 
                   dpi35_vs_Gp_J2_Gp_dpi35_UP, dpi35_vs_Gp_J2_Gp_J2_UP, 
                   dpi35_vs_Gp_MALE_Gp_dpi35_UP, dpi35_vs_Gp_MALE_Gp_MALE_UP, 
                   dpi7_vs_Gp_EGG_Gp_dpi7_UP, dpi7_vs_Gp_EGG_Gp_EGG_UP, 
                   dpi7_vs_Gp_J2_Gp_dpi7_UP, dpi7_vs_Gp_J2_Gp_J2_UP, 
                   dpi7_vs_Gp_MALE_Gp_dpi7_UP, dpi7_vs_Gp_MALE_Gp_MALE_UP, 
                   EGG_vs_Gp_J2_Gp_EGG_UP, EGG_vs_Gp_J2_Gp_J2_UP, 
                   EGG_vs_Gp_MALE_Gp_EGG_UP, EGG_vs_Gp_MALE_Gp_MALE_UP, 
                   J2_vs_Gp_MALE_Gp_J2_UP, J2_vs_Gp_MALE_Gp_MALE_UP, options.out)
