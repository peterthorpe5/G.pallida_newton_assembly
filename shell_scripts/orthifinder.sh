#! -cwd
cd /storage/home/users/pjt6/newton/comparative_genomics/all_nems

conda activate orthofinder

orthofinder -S diamond -t 8 -f /storage/home/users/pjt6/newton/comparative_genomics/all_nems

 python ../mcl_to_cafe.py -sp "GROS Hsc GPLIN GPALN Minc Hetgly" -i ./Results_*/Orthogroups.txt -o nematodes_orthofinder.clusters.txt
 
 
