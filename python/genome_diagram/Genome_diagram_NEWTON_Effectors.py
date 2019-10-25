
"""
this is an untaural process of voodoo to make a very pretty picture!!!

"""

############################################################################
#imports.......

from reportlab.lib import colors
from reportlab.lib.units import cm
# Biopython core
from Bio import SeqIO
from Bio import GenBank
from Bio.SeqFeature import SeqFeature, FeatureLocation
#to loop over files in the folder
import os

# Bio.Graphics.GenomeDiagram
from Bio.Graphics.GenomeDiagram import Diagram

##############################################################################

################################################################

MIN_GAP_JAGGY = 1000
def add_jaggies(contig_seq, offset, gd_contig_features):
    """Add JAGGY features for any run of NNNN or XXXX in sequence."""
    contig_seq = contig_seq.upper().replace("X", "N")
    i = 0
    j = 0
    NNN = "N" * MIN_GAP_JAGGY
    while i < len(contig_seq):
        i = contig_seq.find(NNN, i)
        if i == -1:
            return
        j = i
        while j < len(contig_seq) and contig_seq[j] == "N":
            j += 1
        #print("Adding jaggy")
        gd_contig_features.add_feature(SeqFeature(FeatureLocation(offset+i,
                                                                  offset+j)),
                                       sigil="JAGGY",
                                       color=colors.slategrey,
                                       border=colors.black)
        i = j + 1


def squash_exons(feature):
    """Makes a new SewqFeature discarding the exon information."""
    start = feature.location.start
    end = feature.location.end
    strand = feature.location.strand
    return SeqFeature(FeatureLocation(start, end, strand),
                      type=feature.type,
                      qualifiers=feature.qualifiers)


def draw_me_something_nice (infile, outfile, outfile2):
    """function to draw genome diagrams by looping over
a load of embl files in a folder>>> this is supposed to add
effectors of interest on as coloured items"""
    genbank_entry = SeqIO.read(open(infile), "embl")
    name_for_info_out = infile.split(".embl")[0] + "effecotr_info.txt"
    f_general_output = open(name_for_info_out, "w")
    #print "im here"
    gdd = Diagram('Test Diagram')
    #Add a track of features,
    gdt_features = gdd.new_track(1, greytrack=True,
                                 name="CDS Features",
                                 scale_largetick_interval=100000,
                                 scale_smalltick_interval=5000,
                                 scale_fontsize=3,
                                 scale_format = "SInt",
                                 greytrack_labels=False, #e.g. 5
                                 height=0.75)

    #We'll just use one feature set for these features,
    gds_features = gdt_features.new_set()

    add_jaggies(str(genbank_entry.seq), 0, gds_features)
    
     #genes of interest
    effectors = """#gene
GPALN001111
GPALN001252
GPALN001912
GPALN002106
GPALN002290
GPALN002295
GPALN002300
GPALN002383
GPALN002386
GPALN002387
GPALN002593
GPALN002947
GPALN003010
GPALN003306
GPALN003381
GPALN003415
GPALN003793
GPALN003794
GPALN003795
GPALN003831
GPALN003952
GPALN003970
GPALN003975
GPALN004254
GPALN004470
GPALN004493
GPALN004587
GPALN004712
GPALN004734
GPALN004862
GPALN004897
GPALN005042
GPALN005067
GPALN005090
GPALN005100
GPALN005105
GPALN005801
GPALN005901
GPALN005903
GPALN005905
GPALN005953
GPALN006035
GPALN006057
GPALN006059
GPALN006061
GPALN006067
GPALN006752
GPALN006754
GPALN006755
GPALN006756
GPALN006759
GPALN006766
GPALN006769
GPALN006775
GPALN006818
GPALN006828
GPALN006839
GPALN006853
GPALN006856
GPALN006945
GPALN007181
GPALN007436
GPALN007445
GPALN007670
GPALN007708
GPALN007837
GPALN008101
GPALN009056
GPALN009444
GPALN009458
GPALN009497
GPALN009498
GPALN009532
GPALN009580
GPALN009586
GPALN009589
GPALN009796
GPALN009815
GPALN009825
GPALN009837
GPALN009900
GPALN009918
GPALN010093
GPALN010199
GPALN010232
GPALN010540
GPALN010542
GPALN010554
GPALN010603
GPALN010659
GPALN010702
GPALN010737
GPALN010789
GPALN010793
GPALN010968
GPALN010970
GPALN011399
GPALN011823
GPALN011858
GPALN011865
GPALN012007
GPALN012056
GPALN012064
GPALN012287
GPALN013168
GPALN013277
GPALN013280
GPALN013347
GPALN013348
GPALN013350
GPALN013383
GPALN013384
GPALN013459
GPALN013480
GPALN013496
GPALN014034
GPALN014145
GPALN014146
GPALN014235
GPALN014268
GPALN014324
GPALN014327
GPALN014354
GPALN014355
GPALN014357
GPALN014378
GPALN014379
GPALN014381
GPALN014395
GPALN014477
GPALN014498
GPALN014747
GPALN014750
GPALN014881
GPALN015014
GPALN015073
GPALN015211
GPALN015248
GPALN015280
GPALN015295
GPALN015296
GPALN015298
GPALN015299
GPALN015301
GPALN015302
GPALN015304
GPALN015309
GPALN015425
GPALN015632
GPALN016298
GPALN016343
GPALN016360
GPALN016380
GPALN002204
GPALN002370
GPALN002666
GPALN002991
GPALN003997
GPALN004009
GPALN005554
GPALN007648
GPALN012415
GPALN013387
GPALN013545
GPALN014713
GPALN015272
GPALN015314
GPALN015605
GPALN015654
GPALN002028
GPALN002288
GPALN002969
GPALN003083
GPALN003416
GPALN003852
GPALN004011
GPALN004265
GPALN004480
GPALN004881
GPALN004901
GPALN005017
GPALN005038
GPALN006124
GPALN007079
GPALN007178
GPALN009323
GPALN010231
GPALN010636
GPALN011715
GPALN011857
GPALN012062
GPALN014261
GPALN014271
GPALN014665
GPALN014857
GPALN015013
GPALN015100
GPALN015116
GPALN015172
GPALN015218
GPALN000707
GPALN001149
GPALN001153
GPALN001281
GPALN001284
GPALN001315
GPALN001641
GPALN001729
GPALN001738
GPALN001745
GPALN002294
GPALN002346
GPALN002349
GPALN002377
GPALN002379
GPALN002466
GPALN002494
GPALN003077
GPALN003222
GPALN003326
GPALN003368
GPALN003369
GPALN003846
GPALN003860
GPALN003876
GPALN003882
GPALN003891
GPALN003905
GPALN003908
GPALN003910
GPALN003911
GPALN003912
GPALN003913
GPALN003925
GPALN003942
GPALN003943
GPALN003946
GPALN003949
GPALN003953
GPALN003954
GPALN003955
GPALN003977
GPALN003990
GPALN004007
GPALN004008
GPALN004010
GPALN004014
GPALN004017
GPALN004018
GPALN004064
GPALN004130
GPALN004342
GPALN004369
GPALN004380
GPALN004410
GPALN004411
GPALN004506
GPALN004534
GPALN004553
GPALN004554
GPALN004555
GPALN004557
GPALN004678
GPALN004679
GPALN004681
GPALN004798
GPALN005064
GPALN005129
GPALN005160
GPALN005161
GPALN005611
GPALN005738
GPALN005986
GPALN006031
GPALN006038
GPALN006112
GPALN006223
GPALN006413
GPALN006581
GPALN006596
GPALN006778
GPALN006780
GPALN006782
GPALN006860
GPALN006911
GPALN007051
GPALN007058
GPALN007072
GPALN007120
GPALN007129
GPALN007132
GPALN007139
GPALN007179
GPALN007201
GPALN007443
GPALN007647
GPALN007696
GPALN007748
GPALN007796
GPALN007811
GPALN007848
GPALN007899
GPALN008074
GPALN008097
GPALN008098
GPALN008100
GPALN008102
GPALN008108
GPALN008152
GPALN008161
GPALN008462
GPALN008535
GPALN009441
GPALN009443
GPALN009492
GPALN009505
GPALN009640
GPALN009669
GPALN009670
GPALN009695
GPALN009823
GPALN009839
GPALN009902
GPALN010067
GPALN010126
GPALN010127
GPALN010171
GPALN010316
GPALN010321
GPALN010414
GPALN010416
GPALN010432
GPALN010433
GPALN010511
GPALN010519
GPALN010534
GPALN010536
GPALN010598
GPALN010602
GPALN010621
GPALN010625
GPALN010778
GPALN010795
GPALN010824
GPALN011797
GPALN011812
GPALN011852
GPALN012010
GPALN012025
GPALN012067
GPALN012099
GPALN012127
GPALN012134
GPALN012284
GPALN012357
GPALN012358
GPALN012366
GPALN012838
GPALN013104
GPALN013109
GPALN013144
GPALN013150
GPALN013385
GPALN014005
GPALN014368
GPALN014369
GPALN014370
GPALN014371
GPALN014372
GPALN014377
GPALN014397
GPALN014398
GPALN014514
GPALN014539
GPALN014576
GPALN014672
GPALN014707
GPALN014746
GPALN014851
GPALN014865
GPALN014866
GPALN014867
GPALN014868
GPALN014885
GPALN014904
GPALN015061
GPALN015177
GPALN015178
GPALN015179
GPALN015181
GPALN015182
GPALN015183
GPALN015186
GPALN015188
GPALN015193
GPALN015243
GPALN015262
GPALN015279
GPALN015285
GPALN015291
GPALN015297
GPALN015738
GPALN015769
GPALN016090
GPALN016091
GPALN016117
GPALN016181
GPALN016188
GPALN016330
GPALN016378""".split("\n")
    

    SPRYSEC = """GPALN012056.T1
GPALN009532.T1
GPALN003794.T1
GPALN014357.T1
GPALN010968.T1
GPALN001352.T1
GPALN006035.T1
GPALN007139.T1
GPALN013168.T1
GPALN006853.T1
GPALN010970.T1
GPALN014477.T1
GPALN015302.T1
GPALN012007.T1
GPALN015309.T1
GPALN010793.T1
GPALN006818.T1
GPALN013114.T1
GPALN006860.T1
GPALN009815.T1
GPALN006839.T1
GPALN006856.T1
GPALN004734.T1
GPALN006596.T1
GPALN013383.T1
GPALN011823.T1
GPALN012287.T1
GPALN009918.T1
GPALN014398.T1
GPALN010231.T1
GPALN009669.T1
GPALN010232.T1
GPALN013348.T1
GPALN013350.T1
GPALN010645.T1
GPALN010093.T1
GPALN014397.T1
GPALN002288.T1
GPALN002300.T1
GPALN011858.T1
GPALN015298.T1
GPALN013480.T1
GPALN009458.T1
GPALN010789.T1
GPALN007168.T1
GPALN008646.T1
GPALN006775.T1
GPALN015295.T1
GPALN004897.T1
GPALN002290.T1
GPALN015013.T1
GPALN014271.T1
GPALN015632.T1
GPALN015301.T1
GPALN014355.T1
GPALN007445.T1
GPALN015280.T1
GPALN007711.T1
GPALN015314.T1
GPALN010569.T1
GPALN007132.T1
GPALN006828.T1
GPALN004881.T1
GPALN007129.T1
GPALN013385.T1
GPALN003057.T1
GPALN015407.T1
GPALN004265.T1
GPALN014395.T1
GPALN012062.T1
GPALN001780.T1
GPALN012064.T1
GPALN007120.T1
GPALN005953.T1
GPALN003793.T1
GPALN015813.T1
GPALN016040.T1""".split("\n")
    
    SPRY = """
""".split("\n")
    Dorsal_set = set([])
    J2_set = set([])
    dpi_14_set = set([])
    names = set([])
    effector_list = []
    for i in effectors:
        if i not in names:
            names.add(i + ".T1")
            effector_list.append(i + ".T1")
    dpi_14 = []
    J2 = []
    count = 0
    for feature in genbank_entry.features:
        count = count +1
        shape = "ARROW"
        #if feature.type not in ["CDS", "tRNA", "rRNA"] :
        if feature.type in ["source", "gene"] :#["source", "CDS"]
            #print "CDS"
            #We're going to ignore these (ignore genes as the CDS is enough)
            continue

        #Note that I am using strings for color names, instead
        #of passing in color objects.  This should also work!
        color2 = "grey"
        if feature.type == "tRNA" :
            color = "red"
        elif feature.type == "rRNA":
            color = "purple"
        elif feature.type == "gap":
            color = "grey"
            shape = "JAGGY"
            feature.strand = None #i.e. draw it strandless
        elif feature.type != "CDS" :
            color = "lightgreen"
        # adding two features per gene, so not just odd/even:
        #elif len(gds_features) % 4 == 0 :
        elif count % 2 ==0:
            color = "blue"
            color2 = "lightblue"
            color = colors.Color(0, 0, 1, 0.4)
            color2 = colors.Color(.678431,.847059,.901961,0.2)
        else :
            color = "green"
            color2 = "lightgreen"
            color = colors.Color(0, 0.501961, 0, 0.4)
            color2 = colors.Color(0.564706, 0.933333, 0.564706, 0.2)
        #colour the Dorsal genes yellow
        
        for gene_name in effector_list:
            # print(feature.qualifiers.get("locus_tag", [None])[0].replace(";", ""))
            if feature.qualifiers.get("locus_tag",
                                      [None])[0].replace(";", "") in gene_name.rstrip():
                color = "red"
                color2 = "pink"
                f_general_output.write("effector\t%s\t%s\n" %(infile, gene_name))
                print("effector\t%s\t%s\n" %(infile, gene_name))

        for gene_name in SPRYSEC:
            #print(feature.qualifiers.get("locus_tag",
                                      #[None])[0].replace(";", ""), gene_name)
            if feature.qualifiers.get("locus_tag",
                                      [None])[0].replace(";", "") in gene_name.rstrip():
                color = "blue"
                color2 = "lightblue"
                f_general_output.write("SPRYSEC\t%s\t%s\n" %(infile, gene_name))
                print("SPRYSEC\t%s\t%s\n" %(infile, gene_name))



        gds_features.add_feature(squash_exons(feature), color=color2,
                                 sigil="BOX",
                                 #sigil=shape,
                                 arrowshaft_height=0.8,
                                 arrowhead_length=0.5,
                                 label_position = "start",
                                 label_size = 1,label_angle = 90,
                                 label=True)
        # Don't want the line round the feature as starts to overlap
        gds_features.add_feature(feature, border=False, color=color,
                                 sigil=shape,
                                 arrowshaft_height=0.6,
                                 arrowhead_length=0.5,
                                 label_position = "start",
                                 label_size = 1,label_angle = 90,
                                 label=False)
        #if count/1000.0==3:
            #print count

    

        #And draw it...
    #print "im now drawing it"    
    gdd.draw(format='linear', orientation='landscape',
             tracklines=False, pagesize='A4', fragments=10)
    gdd.write(outfile, 'PDF')
    gdd.write("GROS_linear.svg", 'SVG')

    #And a circular version
    #Change the order and leave an empty space in the center:
    gdd.move_track(1,3)
    gdd.draw(format='circular', tracklines=False, pagesize=(30*cm,30*cm))
    gdd.write(outfile2, 'PDF')
    gdd.write("GROS_circ.svg", 'SVG')

    #print "Done"
    #print "Dorsal \t", Dorsal_set
    #print "J2 \t", J2_set
    #print "dpi_14 \t", dpi_14_set
    #return None #Dorsal_set, J2_set, dpi_14_set

###############################################################################




#to run the code

Dorsal_set2 = set([])
J2_set2 = set([])
dpi_14_set2 = set([])

contigs = """
""".split()

########################################

#load the genbank file that contains the genes
#gbk_filename = "nGr.v1.1.2015_04_09.V2.gbk"
#genbank_entry = SeqIO.read(open(gbk_filename), "genbank")

#genbank_entry = SeqIO.read(open("nGr.v1.1.2015_04_09.gff3.gb"), "genbank")

#print draw_me_something_nice("GROS_00001.embl",\
                             #"GROS_genome_daigram_w_effector20150523.v1.pdf", "out2.out")

for filename in os.listdir("."):
    if not filename.endswith(".embl") : continue
    infile = filename
    outfile = filename.split(".embl")[0] + ".pdf"
    outfile2 = filename.split(".embl")[0] + ".svg"

    draw_me_something_nice (infile, outfile, outfile2)



##for i in contigs:
##    #print i
##    infile = i
##    outfile = "GD_%s.pdf" %(i)
##    #outfile2 = 
##    print draw_me_something_nice (infile, outfile)
##    #print Dorsal_set
##    #Dorsal_set2.add(Dorsal_set)
##    #J2_set2.add(J2_set)
##    #dpi_14_set2.add(dpi_14_set)


    
#f= open ('Dorsal_J2_dpi_14_on_which_contigs.txt','w')
#print >> f, "Dorsal\t", Dorsal_set2
#print >> f, "J2\t", J2_set2
#print >> f, "dpi_14\t", dpi_14_set2

