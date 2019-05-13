#!/bin/bash

cd /data3/Payal/Projects/UDN_Hayk/Patient



hisat2 -t  -p 8 -x /data3/Payal/Genomes/Human/Human_Genome_GRCH38/Hisat2_Human_Primary_Index/Homo_sapiens.GRCh38.dna.primary_assembly \
-1 /data3/Payal/Projects/UDN_Hayk/Patient/FASTQ_Files/Lane1/UDN287643_P_Blood_S1_L001_R1_001.fastq.gz \
-2 /data3/Payal/Projects/UDN_Hayk/Patient/FASTQ_Files/Lane1/UDN287643_P_Blood_S1_L001_R2_001.fastq.gz \
-S /data3/Payal/Projects/UDN_Hayk/Patient/Hisat2_Output/Lane1/UDN287643_P_Blood_S1_L001_001.sam \
--summary-file /data3/Payal/Projects/UDN_Hayk/Patient/Hisat2_Output/Lane1/UDN287643_P_Blood_S1_L001_001_hisat2_aligment_summary.txt

Lib = 1
RGID
Patient Lane1 = 1
Patient Lane2 = 2
Patient Lane3 = 3
Patient Lane4 = 4

Lib = 2
RGID
Father Lane1 = 5
Father Lane2 = 6
Father Lane3 = 7
Father Lane4 = 8

Lib = 3
RGID
Mother Lane1 = 9
Mother Lane2 = 10
Mother Lane3 = 11
Mother Lane4 = 12

picard AddOrReplaceReadGroups I=UDN287643_P_Blood_S1_L001_001.sam O=UDN287643_P_Blood_S1_L001_001_RG.sam RGID=1 RGLB=lib1 RGPL=Illumina RGPU=S1_L001 RGSM=UDN287643_P_Blood

samtools view -H UDN287643_P_Blood_S1_L002_001_RG.sam | grep '@RG'

samtools sort -o UDN287643_P_Blood_S1_L004_001_RG_sorted.bam -O bam UDN287643_P_Blood_S1_L004_001_RG.sam

samtools merge UDN287643_P_Blood_S1.bam /data3/Payal/Projects/UDN_Hayk/Patient/Hisat2_Output/Lane1/UDN287643_P_Blood_S1_L001_001_RG_sorted.bam /data3/Payal/Projects/UDN_Hayk/Patient/Hisat2_Output/Lane2/UDN287643_P_Blood_S1_L002_001_RG_sorted.bam /data3/Payal/Projects/UDN_Hayk/Patient/Hisat2_Output/Lane3/UDN287643_P_Blood_S1_L003_001_RG_sorted.bam /data3/Payal/Projects/UDN_Hayk/Patient/Hisat2_Output/Lane4/UDN287643_P_Blood_S1_L004_001_RG_sorted.bam


rsem-prepare-reference --gtf /data3/Payal/Genomes/Human/Human_Genome_GRCH38/Homo_sapiens.GRCh38.91.gtf --bowtie2 /data3/Payal/Genomes/Human/Human_Genome_GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa /data3/Payal/Genomes/Human/RSEM_Reference_Human/RSEM_Homo_sapiens.GRCh38.dna.primary_assembly

rsem-calculate-expression -p 8 --paired-end --bowtie2 --estimate-rspd --append-names --output-genome-bam \
 /data3/Payal/Projects/UDN_Hayk/Patient/UDN287643_P_Blood_S1_L001_R1_001.fastq.gz,/data3/Payal/Projects/UDN_Hayk/Patient/UDN287643_P_Blood_S1_L002_R1_001.fastq.gz,/data3/Payal/Projects/UDN_Hayk/Patient/UDN287643_P_Blood_S1_L003_R1_001.fastq.gz,/data3/Payal/Projects/UDN_Hayk/Patient/UDN287643_P_Blood_S1_L004_R1_001.fastq.gz \
 /data3/Payal/Projects/UDN_Hayk/Patient/UDN287643_P_Blood_S1_L001_R2_001.fastq.gz,/data3/Payal/Projects/UDN_Hayk/Patient/UDN287643_P_Blood_S1_L002_R2_001.fastq.gz,/data3/Payal/Projects/UDN_Hayk/Patient/UDN287643_P_Blood_S1_L003_R2_001.fastq.gz,/data3/Payal/Projects/UDN_Hayk/Patient/UDN287643_P_Blood_S1_L004_R2_001.fastq.gz\
 /data3/Payal/Genomes/Human/RSEM_Reference_Human/RSEM_Homo_sapiens.GRCh38.dna.primary_assembly /data3/Payal/Projects/UDN_Hayk/Patient/RSEM_Expression_Results/UDN287643_P_Blood_S1_L001

 rsem-calculate-expression -p 8 --paired-end UDN287643_P_Blood_S1_L001_R1_001.fastq,UDN287643_P_Blood_S1_L002_R1_001.fastq,UDN287643_P_Blood_S1_L003_R1_001.fastq,UDN287643_P_Blood_S1_L004_R1_001.fastq UDN287643_P_Blood_S1_L001_R2_001.fastq,UDN287643_P_Blood_S1_L002_R2_001.fastq,UDN287643_P_Blood_S1_L003_R2_001.fastq,UDN287643_P_Blood_S1_L004_R2_001.fastq /data3/Payal/Genomes/Human/RSEM_Reference_Human_GRCh37/RSEM_Homo_sapiens.GRCh37.75.dna.primary_assembly /data3/Payal/Projects/UDN_Hayk/Patient/RSEM_Expression_Results_GRCh37/UDN287643_P
/data3/Payal/Projects/UDN_Hayk/Patient/
 rsem-calculate-expression -p 8 --bowtie2 --paired-end UDN287643_P_Blood_S1_L001_R1_001.fastq,UDN287643_P_Blood_S1_L002_R1_001.fastq,UDN287643_P_Blood_S1_L003_R1_001.fastq,UDN287643_P_Blood_S1_L004_R1_001.fastq UDN287643_P_Blood_S1_L001_R2_001.fastq,UDN287643_P_Blood_S1_L002_R2_001.fastq,UDN287643_P_Blood_S1_L003_R2_001.fastq,UDN287643_P_Blood_S1_L004_R2_001.fastq /data3/Payal/Genomes/Human/RSEM_Reference_Human_GRCh37/RSEM_Homo_sapiens.GRCh37.75.dna.primary_assembly RSEM_Expression_Results_GRCh37/UDN287643_P &