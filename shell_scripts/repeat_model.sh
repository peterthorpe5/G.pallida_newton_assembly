#!/bin/bash
#$ -cwd

cd /storage/home/users/User_name/shelf_apps/newton/newton_wtdgb
conda activate repeatmasking
conda activate python27

python2 run_pipeline.py

