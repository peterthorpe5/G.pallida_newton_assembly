#$ -cwd
cd /storage/home/users/User_name/shelf_apps/newton/newton_wtdgb/L4000
# purger, purged, output all data, blobtools

#python ~/misc_python/convert_file_format/convert_fq_to_fa.py -i raw_reads.fastq -o raw_reads.fasta

conda activate python27

python /shelf/apps/User_name/apps/finishingTool/finisherSC.py -par 8 /storage/home/users/User_name/shelf_apps/newton/newton_wtdgb/L4000/ /shelf/apps/User_name/conda/envs/python27/bin/


