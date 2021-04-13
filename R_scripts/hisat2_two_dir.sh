#!/bin/bash

read1=$(ls /home/pbanerjee/Payal/test_data/read1/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz)
read2=$(ls /home/pbanerjee/Payal/test_data/read2/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz)

ref_dir='/data3/Payal/Genomes/Human/Hisat2_ERCC_Human_Genome_Index/'
ref=$(ls $ref_dir)

#echo "$ref"
#cd "$ref_dir"
#pwd
#ls $ref_dir 

for i in "$ref";
do
#	echo $i
#Homo_sapiens.GRCh38.dna.toplevel_ERCC92.1.ht2l
	name=$(echo $i | cut -c -39)
#	echo $name
	hisat_ref=$ref_dir$name
#	echo $hisat_ref


	echo "The R1 file is $read1";
	echo "The R2 file is  $read2";

	hisat2 --threads 8 -t --summary-file aligned_hisat2.txt -x "$hisat_ref" -1 $read1 -2 $read2 -S aligned_test.sam

done

echo "The output file is aligned_test.sam"       
echo "Finished running Hisat2"
