#!/bin/bash

SAMPLEID=human002

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
$dir_scy/scythe -a adaptors.fa -q sanger -o 1_scy_${SAMPLEID}_pe1.fastq ${SAMPLEID}_pe1.fastq 
$dir_scy/scythe -a adaptors.fa -q sanger -o 1_scy_${SAMPLEID}_pe2.fastq ${SAMPLEID}_pe2.fastq 

##############
### STEP 2 ### SICKLE TO TRIM OFF LOW QUALITY SEQUENCE
$dir_sic/sickle pe -f 1_scy_${SAMPLEID}_pe1.fastq -r 1_scy_${SAMPLEID}_pe2.fastq -t sanger -o 2_sic_${SAMPLEID}_pe1.fastq -p 2_sic_${SAMPLEID}_pe2.fastq -s 2_sic_${SAMPLEID}_sing.fastq -q 20 -l 20

##############
### STEP 3 - BWA map and align create sai file
$dir_bwa/bwa aln $dir_ref/chr14.fa 2_sic_${SAMPLEID}_pe1.fastq > 3_bwa_${SAMPLEID}_pe1.sai 
$dir_bwa/bwa aln $dir_ref/chr14.fa 2_sic_${SAMPLEID}_pe2.fastq > 3_bwa_${SAMPLEID}_pe2.sai 

### STEP 4 - BWA map and align create sam file
$dir_bwa/bwa sampe $dir_ref/chr14.fa 3_bwa_${SAMPLEID}_pe1.sai 3_bwa_${SAMPLEID}_pe2.sai 2_sic_${SAMPLEID}_pe1.fastq 2_sic_${SAMPLEID}_pe2.fastq > 4_bwa_${SAMPLEID}_pe12.sam

#############
### STEP 5 - PICARD TOOLS TO DO FILE FORMATING - SORT BY COORDINATE
java -jar $dir_pic/SortSam.jar I=4_bwa_${SAMPLEID}_pe12.sam O=5_bwasort_${SAMPLEID}_pe12.bam SO=coordinate CREATE_INDEX=true

############
### STEP 6 - SAMTOOLS - FILTER READS THAT ARE UNMAPPED OR HAVE LOW MAPPING QUALITY
$dir_sam/samtools view -F 4 -q 30 -b 5_bwasort_${SAMPLEID}_pe12.bam > 6_bwafil_${SAMPLEID}_pe12.bam

