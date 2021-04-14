# This scipt is adapted by one that Peter C (James Hutton) Wrote.

# Peter Thorpe

from reportlab.lib import colors
from reportlab.lib.units import cm
from reportlab.graphics.shapes import String
from optparse import OptionParser
import sys

from Bio import SeqIO
from Bio.Graphics import BasicChromosome


class LabelledSpacerSegment(BasicChromosome.SpacerSegment):
    def __init__(self, caption=None):
        self.caption = caption
        self.title_size = 12

    def draw(self, cur_drawing):
        if self.caption:
            x_position = 0.5 * (self.end_x_position + self.start_x_position)
            y_position = 0.5 * (self.end_y_position + self.start_y_position)
            label_string = String(x_position, y_position, self.caption)
            label_string.fontName = "Times-Bold"
            label_string.fontSize = self.title_size
            label_string.textAnchor = "middle"
            cur_drawing.add(label_string)


def strip_genes(gene):
    """function to reduce the extar writing around the gene names"""
    gene = gene.split(";")[0]
    gene = gene.replace("ID=", "")
    return gene


def vcf_to_features(vcf_filename, strand, color):
    features = {}
    with open(vcf_filename) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                pass
            if not line.strip():
                continue  # if the last line is blank
            parts = line.split("\t")
            try:
                ref = parts[0]
                pos = int(parts[1])
                end = int(parts[2])
                gene = (parts[3])
                gene = strip_genes(gene)
                strand = "+"
            except ValueError:
                continue
            # (start, end, strand, caption, color, fill_color)
            entry = (pos, end, strand, "", color)
            try:
                features[ref].append(entry)
            except KeyError:
                features[ref] = [entry]
    return features


if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.1.0")
    sys.exit(0)


usage = """Use as follows:
python Check..py -h 

"""

parser = OptionParser(usage=usage)
parser.add_option("-i", "--in",
                  dest="in_fasta",
                  default="Gpal_newton_newton.scaffolds.fa",
                  help="in filename (fasta file)",
                  metavar="FILE")
parser.add_option("-c", "--coordinates",
                  dest="genes_coordinates",
                  default="effectors.coordintaes",
                  help="scaff\tstart\tstop\tGENE\n coodinates file",
                  metavar="FILE")
parser.add_option("-t", "--title",
                  dest="title",
                  default="Scaffold diagrams",
                  help="title that the plots are given",
                  metavar="FILE")
parser.add_option("-o", "--output",
                  dest="out_file",
                  default="effector_families",
                  help="Output filename",
                  metavar="FILE")

(options, args) = parser.parse_args()
outfile = options.out_file
in_fasta = options.in_fasta
genes_coordinates = options.genes_coordinates
title = options.title

if __name__ == '__main__':
    fasta_file = in_fasta
    forward_vcf = genes_coordinates
    forward_color = "red"
    reverse_vcf = "SPRY.coordintaes"
    reverse_color = "white"

    scaffold_count = 16 # 38  # how many columns to draw
    telomere_length = 50000  # For illustration
    spacer_length = 500000
    max_len = 8303752 + 2*telomere_length + spacer_length  # first scaffold

    chr_diagram = BasicChromosome.Organism()
    chr_diagram.page_size = (29.7 * cm, 21 * cm)  # A4 landscape




    forward_features = vcf_to_features("448.coordintaes", +1, "cyan")
    forward2_features = vcf_to_features("GS.coordintaes", +1, "green")
    forward_features.update(forward2_features)
    forward2_features = vcf_to_features("SPRYSEC.coordinates", +1, "red")
    forward_features.update(forward2_features)

    forward2_features = vcf_to_features("1106.coordintaes", +1, "blue")
    forward_features.update(forward2_features)

    forward2_features = vcf_to_features("CLE.coordintaes", +1, "yellow")
    forward_features.update(forward2_features)
    
    reverse_features = vcf_to_features(reverse_vcf, -1, reverse_color)
    scaffolds = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        f_list = forward_features.get(record.id, [])
        r_list = reverse_features.get(record.id, [])
        if not f_list and not r_list:
            continue
        name = record.id.rsplit("_", 1)[1]  # Scaffold number
        scaffolds.append((name, len(record), f_list + r_list))

    pending = scaffolds[:]
    batches = []
    batch = []
    while pending:
        adding = False
        for i, (name, length, features) in enumerate(pending):
            if sum(_[1]+2*telomere_length + spacer_length for _ in batch) + length+2*telomere_length + spacer_length <= max_len:
                # Fits!
                adding = True
                break
        if adding:
            # Add to current column
            batch.append(pending.pop(i))
        else:
            # Full, start new column
            batches.append(batch)
            batch = []
    print(f"Need {len(batches)} columns")

    for index, batch in enumerate(batches):
        name = f"Col {index+1}"
        print(name)
        cur_chromosome = BasicChromosome.Chromosome("")
        # want the same scale used on all, so they can be
        # compared to each other
        cur_chromosome.scale_num = max_len

        for name, length, features in batch:
            start = BasicChromosome.TelomereSegment()
            start.scale = telomere_length
            cur_chromosome.add(start)

            body = BasicChromosome.AnnotatedChromosomeSegment(length, features)
            body.scale = length
            cur_chromosome.add(body)

            end = BasicChromosome.TelomereSegment(inverted=True)
            end.scale = telomere_length
            cur_chromosome.add(end)

            caption = LabelledSpacerSegment(name)
            caption.scale = spacer_length
            cur_chromosome.add(caption)

        # This chromosome is done
        chr_diagram.add(cur_chromosome)

    # Add key
    cur_chromosome = BasicChromosome.Chromosome("")
    cur_chromosome.scale_num = max_len

    start = BasicChromosome.TelomereSegment()
    start.scale = telomere_length
    cur_chromosome.add(start)

    chunk = 1000000 # 1 Mbp
    length = 8303752 # chunk * 8
    features = [(0, 0, None, "0", colors.black, colors.black)]
    for i in range(0,8):
        features.append(
            ((i+1)*chunk, (i+1)*chunk, None, f"{i+1}", colors.black, colors.black)
        )
        block = [i * chunk, (i+1)* chunk, None, ""]
        if i % 2:
            block += [colors.black, colors.lightgrey]
        else:
            block += [colors.black, colors.white]
        features.append(tuple(block))
    body = BasicChromosome.AnnotatedChromosomeSegment(length, features)
    body.scale = length
    cur_chromosome.add(body)

    end = BasicChromosome.TelomereSegment(inverted=True)
    end.scale = telomere_length
    cur_chromosome.add(end)

    caption = LabelledSpacerSegment("Mbp")
    caption.scale = spacer_length
    cur_chromosome.add(caption)

    chr_diagram.add(cur_chromosome)

    chr_diagram.draw(outfile + ".pdf", title)
    chr_diagram.output_format="png"
    chr_diagram.draw(outfile + ".png", title)
