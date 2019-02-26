#!/bin/bash

### SCRIPT FOR CBI_WORKSHOP2

## These are the tool directories for Hoffman2. 
dir_pic="/u/local/apps/picard-tools/current"
dir_sam="/u/local/apps/samtools/0.1.18"
dir_bwa="/u/local/apps/bwa/current"

## we must specify the paths fo these folders
dir_ref=
dir_scy=
dir_sic=

###############
### STEP 1 ### SCYTHE TO TRIM OFF ADAPTORS 
##$dir_scy/scythe -a adaptors.fa -q sanger -o 1_scy_human001_pe1.fastq human001_pe1.fastq 

##############
### STEP 2 ### SICKLE TO TRIM OFF LOW QUALITY SEQUENCE
##$dir_sic/sickle pe -f 1_scy_human001_pe1.fastq -r 1_scy_human001_pe2.fastq -t sanger -o 2_sic_human001_pe1.fastq -p 2_sic_human001_pe2.fastq -s 2_sic_human001_sing.fastq -q 20 -l 20

##############
### STEP 3 - BWA map and align create sai file
##$dir_bwa/bwa aln $dir_ref/chr14.fa 2_sic_human001_pe1.fastq > 3_bwa_human001_pe1.sai 
 
### STEP 4 - BWA map and align create sam file
##$dir_bwa/bwa sampe $dir_ref/chr14.fa 3_bwa_human001_pe1.sai 3_bwa_human001_pe2.sai 2_sic_human001_pe1.fastq 2_sic_human001_pe2.fastq > 4_bwa_human001_pe12.sam

#############
### STEP 5 - PICARD TOOLS TO DO FILE FORMATING - SORT BY COORDINATE

############
### STEP 6 - SAMTOOLS - FILTER READS THAT ARE UNMAPPED OR HAVE LOW MAPPING QUALITY

