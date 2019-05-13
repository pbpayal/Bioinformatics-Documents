#!/bin/bash

#SBATCH -J Patient3
#SBATCH -o Patient3.out
#SBATCH -e Patient3.err
#SBATCH -p short
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -D /lustre/groups/hirokigrp/Ashk/Payal/ABC/Patient3/
#SBATCH --nice=100

module load bwa/0.7.12

bwa mem -M -t 16 /lustre/groups/hirokigrp/Ashk/Payal/ABC/Father/HG37_bwa/Homo_sapiens.GRCh37.dna.toplevel UDN767369-H7WW7ALXX_s5_1_GSLv3-7_26_SL214930.fastq.gz UDN767369-H7WW7ALXX_s5_2_GSLv3-7_26_SL214930.fastq.gz > /lustre/groups/hirokigrp/Ashk/Payal/ABC/Patient3/UDN767369-H7WW7ALXX_s5_GSLv3-7_26_SL214930.sam