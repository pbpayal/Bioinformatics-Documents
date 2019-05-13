#!/bin/bash

mkdir /data3/Payal/Projects/UDN_Hayk/Patient/Hisat2_Output/Lane1

cd /data3/Payal/Projects/UDN_Hayk/Patient/FASTQ_Files/Lane1

for file1 in $(ls *R1_001.fastq.gz)
for file2 in $(ls *R2_001.fastq.gz)
do
	hisat2 -t  -p 8 -x /data3/Payal/Genomes/Human/Human_Genome_GRCH38/Hisat2_Human_Primary_Index/Homo_sapiens.GRCh38.dna.primary_assembly \
	-1 $file1 \
	-2 $file2 \
	-S aligned_$file_R2.sam \
	--summary-file aligned_$file_R2.txt

	mkdir path/folder name
	mv the output files to the new dir


picard AddOrReplaceReadGroups I=UDN287643_P_Blood_S1_L001_001.sam O=UDN287643_P_Blood_S1_L001_001_RG.sam RGID=1 RGLB=lib1 RGPL=Illumina RGPU=S1_L001 RGSM=UDN287643_P_Blood

Patient Lane1 = 1
Patient Lane2 = 2
Patient Lane3 = 3
Patient Lane4 = 4

Father Lane1 = 5
Father Lane2 = 6
Father Lane3 = 7
Father Lane4 = 8

Mother Lane1 = 9
Mother Lane2 = 10
Mother Lane3 = 11
Mother Lane4 = 12

samtools view -H UDN287643_P_Blood_S1_L002_001_RG.sam | grep '@RG'

samtools sort -o UDN287643_P_Blood_S1_L004_001_RG_sorted.bam -O bam UDN287643_P_Blood_S1_L004_001_RG.sam

samtools merge UDN107024_UF_Blood_S2.bam sorted_UDN107024_UF_Blood_S2_L001_001_RG.sam.bam sorted_UDN107024_UF_Blood_S2_L002_001_RG.sam.bam sorted_UDN107024_UF_Blood_S2_L003_001_RG.sam.bam sorted_UDN107024_UF_Blood_S2_L004_001_RG.sam.bam

htseq-count --format bam --additional-attr=gene_name --idattr=gene_id --order pos -t exon UDN287643_P_Blood_S1.bam /data3/Payal/Genomes/Human/Human_Genome_GRCH38/Homo_sapiens.GRCh38.91.gtf > UDN287643_P_Blood_S1_counts.txt



