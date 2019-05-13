#!/bin/bash

#SBATCH -J samtobam_Patient4
#SBATCH -o samtobam_Patient4.out
#SBATCH -e samtobam_Patient4.err
#SBATCH -p short
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -D /lustre/groups/hirokigrp/Ashk/Payal/ABC/Patient4		
#SBATCH --nice=100

module load samtools/1.3.1

samtools view -S -b UDN767369-HWJL2CCXX_s8_GSLv3-7_26_SL214930.sam > UDN767369-HWJL2CCXX_s8_GSLv3-7_26_SL214930.bam