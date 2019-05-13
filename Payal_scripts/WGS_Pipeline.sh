#!/bin/bash
##First the file name has to be pre-processed
##based on parts that would be before and after the read number
for i in $(ls *.fastq.gz | cut -c -11,14- | uniq)  ## Based on File name D3_S1_L001_R1_001.fastq.gz
do
var1=$(echo ${i}| cut -c -10) #Splitting the script to get the first part D3_S1_L001
var2=$(echo ${i}| cut -c 13-) #Splitting the string to get the second part _001.fastq.gz

##aligning using bwa mem. Note we have to change the function used, directory, etc. as and when required.
bwa mem -M -t 16 FastQ/Homo_sapiens.GRCh37.75.dna_rm.primary_assembly.fa $var1'_R1'_$var2 $var1'_R2'_$var2 > $var1'_'$var2'.sam'
done

##Merging the sam files, sorting and indexing to get corresponding merged bam files
samtools merge -r -c 'Final_'$var1'.bam' *.sam
samtools sort -l 5 -O BAM -o 'Final_'$var1'_sort.bam' 'Final_'$var1'.bam'
samtools index 'Final_'$var1'_sort.bam' 'Final_'$var1'_sort.bai'

##Removing Duplicates
java -jar picard.jar MarkDuplicates  REMOVE_DUPLICATES=true I='Final_'$var1'_sort.bam' O='Final_'$var1'_Remooved_Duplicates_sort.bam' M='Final_'$var1'_Remooved_Duplicates_sort.txt'

##Base Recalibaration

##Structural variant and SNP pipelines as applicable

##Printing  the end of the program
echo "end"