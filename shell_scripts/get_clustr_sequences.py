# script to open up a tab MCL output and count the length of the line. The
# amount of seuqnces in each cluster.

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

file_name = 'test.txt'
script_dir = os.path.dirname(os.path.abspath(__file__))
dest_dir = os.path.join(script_dir, 'effectors')
try:
    os.makedirs(dest_dir)
except OSError:
    print("already exists")
dest_dir = os.path.join(script_dir, 'sub')
try:
    os.makedirs(dest_dir)
except OSError:
    print("already exists")
dest_dir = os.path.join(script_dir, 'busco')
try:
    os.makedirs(dest_dir)
except OSError:
    print("already exists")


def parse_tab_file_get_clusters(fasta_filename,
                                cluster_file, sub_set,
                                effector_set,
                                busco_set,
                                gros_effector):
    """#script to open up a tab MCL output and count the length of the line. The
    # amount of seuqnces in each cluster"""

    MCL_file = open (cluster_file, "r")
    predicted_protein =  SeqIO.index(fasta_filename, "fasta")
    count = int(0)
    for line in MCL_file:
        effector_cluster = ""
        sub_cluster = ""
        busco_cluster = ""
        gros_cluster = ""
        effector_found = False
        sub_found = False
        busco_found = False
        gros_found = False

        MCL_cluster_line = line.rstrip("\n").split()
        cluster_size = len(MCL_cluster_line)
        count += 1
        filename_effectors = "./effectors/effectors_cluster_%d_len_%d.fasta" %(count,cluster_size)
        filename_subs = "./sub/sub_cluster_%d_len_%d.fasta" %(count,cluster_size)
        filename_buscos = "./busco/busco_cluster_%d_len_%d.fasta" %(count,cluster_size)
        filename_gros = "./gros/gros_cluster_%d_len_%d.fasta" %(count,cluster_size)

        nameset = set()
        for i in MCL_cluster_line:
            if i.startswith("OG"):
               continue
            if not len(MCL_cluster_line) > 3:
                continue
            if i in effector_set:
                print("Effector\t%s\t%d" % (i, count))
                effector_found = True
            if i in sub_set:
                print("sub\t%s\t%d" % (i, count))
                sub_found= True
            if i in busco_set:
                print("busco\t%s\t%d" % (i, count))
                busco_found = True
            if i in gros_effector:
                print("GROS\t%s\t%d" % (i, count))
                gros_found = True
        nameset = set()        
        if effector_found:
            for entry in MCL_cluster_line:
                if entry.startswith("OG"):
                    continue
                seq_record = predicted_protein[entry]
                seq_record.description = ""
                if "GROS" in seq_record.id:
                    seq_record.id = seq_record.id.split(".t")[0]
                if seq_record.id not in nameset:
                    info = ">%s\n%s\n" %(seq_record.id, seq_record.seq)
                    effector_cluster = effector_cluster + info
                    nameset.add(seq_record.id)
            filename_eff = open(filename_effectors, "w")
            filename_eff.write(effector_cluster)
            filename_eff.close()
        nameset = set()        
        if sub_found:
            for entry in MCL_cluster_line:
                if entry.startswith("OG"):
                    continue
                seq_record = predicted_protein[entry]
                seq_record.description = ""
                if "GROS" in seq_record.id:
                    seq_record.id = seq_record.id.split(".t")[0]
                if seq_record.id not in nameset:
                    info = ">%s\n%s\n" %(seq_record.id, seq_record.seq)
                    nameset.add(seq_record.id)
                    sub_cluster = sub_cluster + info
            filename_su = open(filename_subs, "w")
            filename_su.write(sub_cluster)
            filename_su.close()
        nameset = set()
        if busco_found:
            for entry in MCL_cluster_line:
                if entry.startswith("OG"):
                    continue
                seq_record = predicted_protein[entry]
                seq_record.description = ""
                if "GROS" in seq_record.id:
                    seq_record.id = seq_record.id.split(".t")[0]
                if seq_record.id not in nameset:
                    info = ">%s\n%s\n" %(seq_record.id, seq_record.seq)
                    nameset.add(seq_record.id)
                    busco_cluster = busco_cluster + info
            filename_b = open(filename_buscos, "w")
            filename_b.write(busco_cluster)
            filename_b.close()
        nameset = set()
        if gros_found:
            for entry in MCL_cluster_line:
                if entry.startswith("OG"):
                    continue
                seq_record = predicted_protein[entry]
                seq_record.description = ""
                if "GROS" in seq_record.id:
                    seq_record.id = seq_record.id.split(".t")[0]
                if seq_record.id not in nameset:
                    info = ">%s\n%s\n" %(seq_record.id, seq_record.seq)
                    nameset.add(seq_record.id)
                    gros_cluster = gros_cluster + info
            filename_g = open(filename_gros, "w")
            filename_g.write(gros_cluster)
            filename_g.close()


            
effectors = set("""#gene
GPALN_001111
GPALN_001252
GPALN_001912
GPALN_002106
GPALN_002290
GPALN_002295
GPALN_002300
GPALN_002383
GPALN_002386
GPALN_002387
GPALN_002593
GPALN_002947
GPALN_003010
GPALN_003306
GPALN_003381
GPALN_003415
GPALN_003793
GPALN_003794
GPALN_003795
GPALN_003831
GPALN_003952
GPALN_003970
GPALN_003975
GPALN_004254
GPALN_004470
GPALN_004493
GPALN_004587
GPALN_004712
GPALN_004734
GPALN_004862
GPALN_004897
GPALN_005042
GPALN_005067
GPALN_005090
GPALN_005100
GPALN_005105
GPALN_005801
GPALN_005901
GPALN_005903
GPALN_005905
GPALN_005953
GPALN_006035
GPALN_006057
GPALN_006059
GPALN_006061
GPALN_006067
GPALN_006752
GPALN_006754
GPALN_006755
GPALN_006756
GPALN_006759
GPALN_006766
GPALN_006769
GPALN_006775
GPALN_006818
GPALN_006828
GPALN_006839
GPALN_006853
GPALN_006856
GPALN_006945
GPALN_007181
GPALN_007436
GPALN_007445
GPALN_007670
GPALN_007708
GPALN_007837
GPALN_008101
GPALN_009056
GPALN_009444
GPALN_009458
GPALN_009497
GPALN_009498
GPALN_009532
GPALN_009580
GPALN_009586
GPALN_009589
GPALN_009796
GPALN_009815
GPALN_009825
GPALN_009837
GPALN_009900
GPALN_009918
GPALN_010093
GPALN_010199
GPALN_010232
GPALN_010540
GPALN_010542
GPALN_010554
GPALN_010603
GPALN_010659
GPALN_010702
GPALN_010737
GPALN_010789
GPALN_010793
GPALN_010968
GPALN_010970
GPALN_011399
GPALN_011823
GPALN_011858
GPALN_011865
GPALN_012007
GPALN_012056
GPALN_012064
GPALN_012287
GPALN_013168
GPALN_013277
GPALN_013280
GPALN_013347
GPALN_013348
GPALN_013350
GPALN_013383
GPALN_013384
GPALN_013459
GPALN_013480
GPALN_013496
GPALN_014034
GPALN_014145
GPALN_014146
GPALN_014235
GPALN_014268
GPALN_014324
GPALN_014327
GPALN_014354
GPALN_014355
GPALN_014357
GPALN_014378
GPALN_014379
GPALN_014381
GPALN_014395
GPALN_014477
GPALN_014498
GPALN_014747
GPALN_014750
GPALN_014881
GPALN_015014
GPALN_015073
GPALN_015211
GPALN_015248
GPALN_015280
GPALN_015295
GPALN_015296
GPALN_015298
GPALN_015299
GPALN_015301
GPALN_015302
GPALN_015304
GPALN_015309
GPALN_015425
GPALN_015632
GPALN_016298
GPALN_016343
GPALN_016360
GPALN_016380
GPALN_002204
GPALN_002370
GPALN_002666
GPALN_002991
GPALN_003997
GPALN_004009
GPALN_005554
GPALN_007648
GPALN_012415
GPALN_013387
GPALN_013545
GPALN_014713
GPALN_015272
GPALN_015314
GPALN_015605
GPALN_015654
GPALN_002028
GPALN_002288
GPALN_002969
GPALN_003083
GPALN_003416
GPALN_003852
GPALN_004011
GPALN_004265
GPALN_004480
GPALN_004881
GPALN_004901
GPALN_005017
GPALN_005038
GPALN_006124
GPALN_007079
GPALN_007178
GPALN_009323
GPALN_010231
GPALN_010636
GPALN_011715
GPALN_011857
GPALN_012062
GPALN_014261
GPALN_014271
GPALN_014665
GPALN_014857
GPALN_015013
GPALN_015100
GPALN_015116
GPALN_015172
GPALN_015218
GPALN_000707
GPALN_001149
GPALN_001153
GPALN_001281
GPALN_001284
GPALN_001315
GPALN_001641
GPALN_001729
GPALN_001738
GPALN_001745
GPALN_002294
GPALN_002346
GPALN_002349
GPALN_002377
GPALN_002379
GPALN_002466
GPALN_002494
GPALN_003077
GPALN_003222
GPALN_003326
GPALN_003368
GPALN_003369
GPALN_003846
GPALN_003860
GPALN_003876
GPALN_003882
GPALN_003891
GPALN_003905
GPALN_003908
GPALN_003910
GPALN_003911
GPALN_003912
GPALN_003913
GPALN_003925
GPALN_003942
GPALN_003943
GPALN_003946
GPALN_003949
GPALN_003953
GPALN_003954
GPALN_003955
GPALN_003977
GPALN_003990
GPALN_004007
GPALN_004008
GPALN_004010
GPALN_004014
GPALN_004017
GPALN_004018
GPALN_004064
GPALN_004130
GPALN_004342
GPALN_004369
GPALN_004380
GPALN_004410
GPALN_004411
GPALN_004506
GPALN_004534
GPALN_004553
GPALN_004554
GPALN_004555
GPALN_004557
GPALN_004678
GPALN_004679
GPALN_004681
GPALN_004798
GPALN_005064
GPALN_005129
GPALN_005160
GPALN_005161
GPALN_005611
GPALN_005738
GPALN_005986
GPALN_006031
GPALN_006038
GPALN_006112
GPALN_006223
GPALN_006413
GPALN_006581
GPALN_006596
GPALN_006778
GPALN_006780
GPALN_006782
GPALN_006860
GPALN_006911
GPALN_007051
GPALN_007058
GPALN_007072
GPALN_007120
GPALN_007129
GPALN_007132
GPALN_007139
GPALN_007179
GPALN_007201
GPALN_007443
GPALN_007647
GPALN_007696
GPALN_007748
GPALN_007796
GPALN_007811
GPALN_007848
GPALN_007899
GPALN_008074
GPALN_008097
GPALN_008098
GPALN_008100
GPALN_008102
GPALN_008108
GPALN_008152
GPALN_008161
GPALN_008462
GPALN_008535
GPALN_009441
GPALN_009443
GPALN_009492
GPALN_009505
GPALN_009640
GPALN_009669
GPALN_009670
GPALN_009695
GPALN_009823
GPALN_009839
GPALN_009902
GPALN_010067
GPALN_010126
GPALN_010127
GPALN_010171
GPALN_010316
GPALN_010321
GPALN_010414
GPALN_010416
GPALN_010432
GPALN_010433
GPALN_010511
GPALN_010519
GPALN_010534
GPALN_010536
GPALN_010598
GPALN_010602
GPALN_010621
GPALN_010625
GPALN_010778
GPALN_010795
GPALN_010824
GPALN_011797
GPALN_011812
GPALN_011852
GPALN_012010
GPALN_012025
GPALN_012067
GPALN_012099
GPALN_012127
GPALN_012134
GPALN_012284
GPALN_012357
GPALN_012358
GPALN_012366
GPALN_012838
GPALN_013104
GPALN_013109
GPALN_013144
GPALN_013150
GPALN_013385
GPALN_014005
GPALN_014368
GPALN_014369
GPALN_014370
GPALN_014371
GPALN_014372
GPALN_014377
GPALN_014397
GPALN_014398
GPALN_014514
GPALN_014539
GPALN_014576
GPALN_014672
GPALN_014707
GPALN_014746
GPALN_014851
GPALN_014865
GPALN_014866
GPALN_014867
GPALN_014868
GPALN_014885
GPALN_014904
GPALN_015061
GPALN_015177
GPALN_015178
GPALN_015179
GPALN_015181
GPALN_015182
GPALN_015183
GPALN_015186
GPALN_015188
GPALN_015193
GPALN_015243
GPALN_015262
GPALN_015279
GPALN_015285
GPALN_015291
GPALN_015297
GPALN_015738
GPALN_015769
GPALN_016090
GPALN_016091
GPALN_016117
GPALN_016181
GPALN_016188
GPALN_016330
GPALN_016378""".split())

busco_set = set("""GPALN_009154-T1
GPALN_001471-T1
GPALN_013430-T1
GPALN_006219-T1
GPALN_012789-T1
GPALN_011026-T1
GPALN_002100-T1
GPALN_011053-T1
GPALN_000825-T1
GPALN_001998-T1
GPALN_011181-T1
GPALN_012686-T1
GPALN_003602-T1
GPALN_005329-T1
GPALN_008212-T1
GPALN_008197-T1
GPALN_006186-T1
GPALN_015726-T1
GPALN_003557-T1
GPALN_014450-T1
GPALN_005127-T1
GPALN_008715-T1
GPALN_008715-T2
GPALN_002803-T1
GPALN_006148-T1
GPALN_005294-T1
GPALN_000523-T1
GPALN_000836-T1
GPALN_003556-T1
GPALN_000444-T1
GPALN_001526-T1
GPALN_015544-T1
GPALN_005318-T1
GPALN_001073-T1
GPALN_008968-T1
GPALN_009531-T1
GPALN_013978-T1
GPALN_001927-T1
GPALN_012984-T1
GPALN_008733-T1
GPALN_005653-T1
GPALN_008825-T1
GPALN_001988-T1
GPALN_014432-T1
GPALN_011477-T1
GPALN_011128-T1
GPALN_006192-T1
GPALN_012747-T1
GPALN_013740-T1
GPALN_016206-T1
GPALN_015032-T1
GPALN_005684-T1
GPALN_011422-T1
GPALN_000639-T1
GPALN_010786-T1
GPALN_010787-T1
GPALN_010792-T1
GPALN_006499-T1
GPALN_000413-T1
GPALN_013870-T1
GPALN_002851-T1
GPALN_005222-T1
GPALN_015998-T1
GPALN_012605-T1
GPALN_012605-T2
GPALN_014737-T1
GPALN_009295-T1
GPALN_009295-T2
GPALN_002937-T1
GPALN_008923-T1
GPALN_000534-T1
GPALN_003489-T1
GPALN_003041-T1
GPALN_003721-T1
GPALN_003721-T2
GPALN_003721-T3
GPALN_008890-T1
GPALN_008898-T1
GPALN_009009-T1
GPALN_002800-T1
GPALN_000057-T1
GPALN_001483-T1
GPALN_008704-T1
GPALN_002618-T1
GPALN_000592-T1
GPALN_011262-T1
GPALN_012519-T1
GPALN_010014-T1
GPALN_002932-T1
GPALN_015719-T1
GPALN_004966-T1
GPALN_012733-T1
GPALN_011986-T1
GPALN_009090-T1
GPALN_000356-T1
GPALN_001400-T1
GPALN_003686-T1
GPALN_010081-T1
GPALN_005875-T1
GPALN_015549-T1
GPALN_000856-T1
GPALN_003661-T1
GPALN_000154-T1
GPALN_008921-T1
GPALN_015459-T1
GPALN_005258-T1
GPALN_000804-T1
GPALN_008911-T1
GPALN_001404-T1
GPALN_009423-T1
GPALN_006604-T1
GPALN_008190-T1
GPALN_011206-T1
GPALN_003685-T1
GPALN_011226-T1
GPALN_011226-T2
GPALN_006341-T1
GPALN_006341-T2
GPALN_003747-T1
GPALN_011367-T1
GPALN_011283-T1
GPALN_012987-T1
GPALN_013742-T1
GPALN_005628-T1
GPALN_007240-T1
GPALN_008750-T1
GPALN_003588-T1
GPALN_014712-T1
GPALN_012745-T1
GPALN_014628-T1
GPALN_013826-T1
GPALN_006271-T1
GPALN_008835-T1
GPALN_000950-T1
GPALN_002825-T1
GPALN_001033-T1
GPALN_002906-T1
GPALN_005438-T1
GPALN_013855-T1
GPALN_003089-T1
GPALN_010761-T1
GPALN_005539-T1
GPALN_001428-T1
GPALN_011251-T1
GPALN_015765-T1
GPALN_012805-T1
GPALN_011356-T1
GPALN_013775-T1
GPALN_007245-T1
GPALN_008851-T1
GPALN_006530-T1
GPALN_012729-T1
GPALN_012994-T1
GPALN_016295-T1
GPALN_005316-T1
GPALN_005285-T1
GPALN_001443-T1
GPALN_001865-T1
GPALN_003570-T1
GPALN_014464-T1
GPALN_012902-T1
GPALN_004994-T1
GPALN_012584-T1
GPALN_009959-T1
GPALN_000072-T1
GPALN_012695-T1
GPALN_003490-T1
GPALN_011052-T1
GPALN_012385-T1
GPALN_014646-T1
GPALN_014791-T1
GPALN_012849-T1
GPALN_008965-T1
GPALN_001932-T1
GPALN_003638-T1
GPALN_013650-T1
GPALN_013759-T1
GPALN_005673-T1
GPALN_012697-T1
GPALN_012700-T1
GPALN_011047-T1
GPALN_010762-T1
GPALN_006336-T1
GPALN_003543-T1
GPALN_001695-T1
GPALN_007949-T1
GPALN_003488-T1
GPALN_003530-T1
GPALN_012953-T1
GPALN_015744-T1
GPALN_001537-T1
GPALN_013699-T1
GPALN_001392-T1
GPALN_015620-T1
GPALN_008428-T1
GPALN_011107-T1
GPALN_003453-T1
GPALN_002566-T1
GPALN_001474-T1
GPALN_002745-T1
GPALN_011027-T1
GPALN_000848-T1
GPALN_005503-T1
GPALN_013790-T1
GPALN_015779-T1
GPALN_003154-T1
GPALN_010974-T1
GPALN_015946-T1
GPALN_002923-T1
GPALN_006448-T1
GPALN_006463-T1
GPALN_013062-T1
GPALN_000875-T1
GPALN_009697-T1
GPALN_000905-T1
GPALN_000295-T1
GPALN_012772-T1
GPALN_000611-T1
GPALN_011469-T1
GPALN_001375-T1
GPALN_008196-T1
GPALN_002786-T1
GPALN_013902-T1
GPALN_014807-T1
GPALN_008765-T1
GPALN_005523-T1
GPALN_011168-T1
GPALN_011273-T1
GPALN_000357-T1
GPALN_005403-T1
GPALN_010862-T1
GPALN_001633-T1
GPALN_015490-T1
GPALN_001587-T1
GPALN_013036-T1
GPALN_002389-T1
GPALN_016081-T1
GPALN_008418-T1
GPALN_011974-T1
GPALN_013953-T1
GPALN_008422-T1
GPALN_010849-T1
GPALN_013679-T1
GPALN_000256-T1
GPALN_011049-T1
GPALN_011666-T1
GPALN_013730-T1
GPALN_011144-T1
GPALN_013038-T1
GPALN_013061-T1
GPALN_015475-T1
GPALN_012781-T1
GPALN_003495-T1
GPALN_007478-T1
GPALN_006180-T1
GPALN_001324-T1
GPALN_009345-T1
GPALN_009192-T1
GPALN_010854-T1
GPALN_005687-T1
GPALN_010941-T1
GPALN_003620-T1
GPALN_001599-T1
GPALN_001667-T1
GPALN_008897-T1
GPALN_016103-T1
GPALN_001566-T1
GPALN_012645-T1
GPALN_003569-T1
GPALN_008709-T1
GPALN_008712-T1
GPALN_015530-T1
GPALN_012904-T1
GPALN_011340-T1
GPALN_008416-T1
GPALN_003573-T1
GPALN_016009-T1
GPALN_008395-T1
GPALN_016011-T1
GPALN_001861-T1
GPALN_011151-T1
GPALN_000417-T1
GPALN_000420-T1
GPALN_013726-T1
GPALN_008621-T1
GPALN_000249-T1
GPALN_002414-T1
GPALN_005331-T1
GPALN_012897-T1
GPALN_006134-T1
GPALN_006309-T1
GPALN_009220-T1
GPALN_001605-T1
GPALN_008768-T1
GPALN_000907-T1
GPALN_000650-T1
GPALN_010004-T1
GPALN_001597-T1
GPALN_013640-T1
GPALN_015985-T1
GPALN_003624-T1
GPALN_010032-T1
GPALN_000479-T1
GPALN_001670-T1
GPALN_000910-T1
GPALN_001661-T1
GPALN_007486-T1
GPALN_015998-T1
GPALN_012460-T1
GPALN_001704-T1
GPALN_003525-T1
GPALN_013894-T1
GPALN_016057-T1
GPALN_000772-T1
GPALN_008829-T1
GPALN_008713-T1
GPALN_000268-T1
GPALN_000968-T1
GPALN_012951-T1
GPALN_001886-T1
GPALN_014721-T1""".split())

sub = set("""GPALN_001169
GPALN_001206
GPALN_001525
GPALN_002386
GPALN_002433
GPALN_002593
GPALN_008104
GPALN_003797
GPALN_003941
GPALN_004056
GPALN_004293
GPALN_004294
GPALN_004460
GPALN_004487
GPALN_004592
GPALN_004595
GPALN_004869
GPALN_005081
GPALN_005082
GPALN_005088
GPALN_005090
GPALN_005097
GPALN_005100
GPALN_005101
GPALN_005105
GPALN_005107
GPALN_005809
GPALN_005901
GPALN_005903
GPALN_005904
GPALN_006022
GPALN_006023
GPALN_006026
GPALN_006099
GPALN_006304
GPALN_006652
GPALN_006654
GPALN_006949
GPALN_007176
GPALN_007507
GPALN_008079
GPALN_008104
GPALN_008135
GPALN_008137
GPALN_008860
GPALN_008886
GPALN_009628
GPALN_009649
GPALN_010131
GPALN_010257
GPALN_010300
GPALN_010417
GPALN_010603
GPALN_010702
GPALN_012301
GPALN_012703
GPALN_013039
GPALN_013946
GPALN_014068
GPALN_014324
GPALN_015174
GPALN_015222
GPALN_015288
GPALN_015365
GPALN_015697
GPALN_015826
GPALN_015875
GPALN_016116
GPALN_016343""".split())


gros = set("""sseqid
GROS_g00017.t1
GROS_g00601.t1
GROS_g00774.t1
GROS_g01784.t1
GROS_g02054.t1
GROS_g02107.t1
GROS_g02358.t1
GROS_g02394.t1
GROS_g02441.t1
GROS_g02469.t1
GROS_g02470.t1
GROS_g02635.t1
GROS_g02798.t1
GROS_g03169.t1
GROS_g03475.t1
GROS_g03476.t1
GROS_g03975.t1
GROS_g04216.t1
GROS_g04366.t1
GROS_g04620.t1
GROS_g04623.t1
GROS_g04662.t1
GROS_g04677.t1
GROS_g04939.t1
GROS_g05241.t1
GROS_g05352.t1
GROS_g05354.t1
GROS_g05390.t1
GROS_g05398.t1
GROS_g05682.t1
GROS_g05707.t1
GROS_g05724.t1
GROS_g05961.t1
GROS_g05962.t1
GROS_g05968.t1
GROS_g05982.t1
GROS_g05988.t1
GROS_g06034.t1
GROS_g06220.t1
GROS_g06322.t1
GROS_g06357.t1
GROS_g06362.t1
GROS_g06363.t1
GROS_g06364.t1
GROS_g06661.t1
GROS_g06952.t1
GROS_g07338.t1
GROS_g07341.t1
GROS_g07372.t1
GROS_g07463.t1
GROS_g07677.t1
GROS_g07949.t1
GROS_g07968.t1
GROS_g08030.t1
GROS_g08139.t1
GROS_g08163.t1
GROS_g08166.t1
GROS_g08169.t1
GROS_g08190.t1
GROS_g08588.t1
GROS_g08634.t1
GROS_g08694.t1
GROS_g08876.t1
GROS_g08879.t1
GROS_g08893.t1
GROS_g08949.t1
GROS_g08953.t1
GROS_g08990.t1
GROS_g09018.t1
GROS_g09038.t1
GROS_g09513.t1
GROS_g09514.t1
GROS_g09568.t1
GROS_g09598.t1
GROS_g09735.t1
GROS_g09961.t1
GROS_g10304.t1
GROS_g10395.t1
GROS_g10505.t1
GROS_g10585.t1
GROS_g10725.t1
GROS_g10726.t1
GROS_g10807.t1
GROS_g10874.t1
GROS_g10924.t1
GROS_g10954.t1
GROS_g11040.t1
GROS_g11397.t1
GROS_g11726.t1
GROS_g11727.t1
GROS_g11775.t1
GROS_g11776.t1
GROS_g11812.t1
GROS_g11888.t1
GROS_g11889.t1
GROS_g12027.t1
GROS_g12028.t1
GROS_g12030.t1
GROS_g12196.t1
GROS_g12320.t1
GROS_g12422.t1
GROS_g12477.t1
GROS_g12501.t1
GROS_g12570.t1
GROS_g12651.t1
GROS_g12705.t1
GROS_g12709.t1
GROS_g12791.t1
GROS_g12817.t1
GROS_g12818.t1
GROS_g12827.t1
GROS_g12966.t1
GROS_g13055.t1
GROS_g13195.t1
GROS_g13374.t1
GROS_g13394.t1
GROS_g13503.t1
GROS_g13646.t1
GROS_g13797.t1
GROS_g13840.t1
GROS_g14123.t1
GROS_g14125.t1
GROS_g14130.t1
GROS_g14149.t1
GROS_g14157.t1
GROS_g14158.t1
GROS_g14167.t1
GROS_g14168.t1
GROS_g14178.t1
GROS_g14179.t1
GROS_g14180.t1
GROS_g14189.t1
GROS_g14194.t1
GROS_g14220.t1
GROS_g14228.t1
GROS_g14232.t1
GROS_g14234.t1
GROS_g14236.t1
GROS_g14275.t1
GROS_g14299.t1
GROS_g14300.t1
GROS_g14306.t1
GROS_g14309.t1
""".split())


effector_set = set([])
for i in effectors:
    i = i + "-T1"
    effector_set.add(i)


sub_set = set([])
for i in sub:
    i = i + "-T1"
    sub_set.add(i)


gros_effector = set([])
for i in gros:
    i = i.split(".t")[0]
    gros_effector.add(i)

    
# busco_set
parse_tab_file_get_clusters("/storage/home/users/pjt6/newton/comparative_genomics/all_nems/all_prot.fasta",
                            "/storage/home/users/pjt6/newton/comparative_genomics/all_nems/Results_Aug12/Orthogroups.txt",
                            sub_set, effector_set, busco_set, gros_effector)





                            
