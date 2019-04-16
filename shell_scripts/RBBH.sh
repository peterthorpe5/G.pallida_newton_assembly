#!/bin/bash
#$ -cwd

cd /storage/home/users/User_name/newton/RBBH


python Blast_RBH_two_fasta_file_evalue.py -a prot --threads 8 -o published_vs_newton_RBBH.tab  Gpal_newton_newton.proteins.fa Gpal.v1.0.AA.fa



