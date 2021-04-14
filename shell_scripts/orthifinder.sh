#! -cwd
cd /storage/home/users/pjt6/newton/comparative_genomics/all_nems

conda activate orthofinder

orthofinder -S diamond -t 8 -f /storage/home/users/pjt6/newton/comparative_genomics/all_nems

 python ../mcl_to_cafe.py -sp "GROS Hsc GPLIN GPALN Minc Hetgly" -i ./Results_*/Orthogroups.txt -o nematodes_orthofinder.clusters.txt
 
 
orthofinder -S diamond -t 32 -f /storage/home/users/John_Jones/newton/effector_finding/

python /storage/home/users/John_Jones/newton/effector_finding/mcl_to_cafe.py -sp "GPALN GPLIN GH5 GH43 GH53 CBM Expansin " -i ./Results_*/Orthogroups.txt -o GPLA_N_eff_HGT_clusters_orthofinder.clusters.txt




python /storage/home/users/John_Jones/newton/effector_finding/mcl_to_cafe.py -sp "GPALN GPLIN GH5 GH43 GH53 CBM Expansin Hgly_G18H08 Chorismate_mutase Hgly_G19B10 Hgly_G4G05 Hgsec11 G._pallida_IA7 RKN_GCprotien28 Hgly_33A09 Hgly_30G12 448 Hgly_G17G01 30G12 E9 Hgly_G23G11 A42 Hgsec3 G8A07_ Hgsec12 IA7 Hgg17 hgsec8 dgl1 IVG9 scn1120. G16H02 H.avenae C52_protein-like Hgly_G12H04 AY135365 Hgly_G20E03 10C02 GpUBI-EP Hgg-20 G10A06 1106 G19C07 747 Hgly_G16A01 G7E05 CLE Hgly_Hgg_20 Invertase CM"  -i ./Results_*/Orthogroups.txt -o EFF_names_as_clusters_HGT_clusters_orthofinder.clusters.txt


python /storage/home/users/John_Jones/newton/effector_finding/mcl_to_cafe.py -sp "GPALN GPLIN GH5 GH43 GH53 CBM Expansin Hgly_G18H08 Chorismate_mutase Hgly_G19B10 Hgly_G4G05 Hgsec11 G._pallida_IA7 RKN_GCprotien28 Hgly_33A09 Hgly_30G12 448 Hgly_G17G01 30G12 E9 Hgly_G23G11 A42 Hgsec3 G8A07_ Hgsec12 IA7 Hgg17 hgsec8 dgl1 IVG9 scn1120. G16H02 H.avenae C52_protein-like Hgly_G12H04 AY135365 Hgly_G20E03 10C02 GpUBI-EP Hgg-20 G10A06 1106 G19C07 747 Hgly_G16A01 G7E05 CLE Hgly_Hgg_20 Invertase CM"  -i Orthogroups_secreted.txt -o EFF_names_as_clusters_HGT_clusters_orthofinder.NEWT_sectred.clusters.txt

