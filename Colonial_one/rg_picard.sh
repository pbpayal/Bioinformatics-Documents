#!/bin/bash

#SBATCH -J rg_Father
#SBATCH -o rg_Father.out
#SBATCH -e rg_Father.err
#SBATCH -p short
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -D /lustre/groups/hirokigrp/Ashk/Payal/ABC/Father	
#SBATCH --nice=100

module load picard-tools/2.6.0

java -jar picard.jar AddOrReplaceReadGroups \
      I=UDN413772-H7WW7ALXX_s5_GSLv3-7_27_SL214931.bam \
      O=UDN413772-H7WW7ALXX_s5_GSLv3-7_27_SL214931_sorted_RG.bam \
      RGID=1 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=UDN413772-H7WW7ALXX \
      RGSM=UDN413772










