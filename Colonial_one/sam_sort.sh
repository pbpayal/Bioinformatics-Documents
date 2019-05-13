#!/bin/bash

#SBATCH -J samsort_Father
#SBATCH -o samtobam_Father.out
#SBATCH -e samtobam_Father.err
#SBATCH -p 128gb
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -D /lustre/groups/hirokigrp/Ashk/Payal/ABC/Father	
#SBATCH --nice=100

module load samtools/1.3.1

for file in $(ls *bam)
do
samtools sort -o sorted$file -O bam $file
done
