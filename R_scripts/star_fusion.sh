#!/bin/bash

#SBATCH -J STAR_Fusion
#SBATCH -o STAR_Fusion.out
#SBATCH -e STAR_Fusion.err
#SBATCH -p short
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -D /lustre/groups/hirokigrp/Ashk/Payal/TS20/STAR_FUSION
#SBATCH --nice=100

module load star/2.5

STAR-Fusion --chimeric_out_sam /lustre/groups/hirokigrp/Ashk/Payal/TS20/Chimeric.out.sam --chimeric_junction /lustre/groups/hirokigrp/Ashk/Payal/TS20/Chimeric.out.junction --genomeDir /lustre/groups/hirokigrp/Ashk/Payal/Human_ERCC_Genome/  --output_dir /lustre/groups/hirokigrp/Ashk/Payal/TS20/STAR_FUSION
