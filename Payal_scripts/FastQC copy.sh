#!/bin/bash
##First the file name has to be pre-processed
##based on parts that would be before and after the read number
for i in $(ls *.fastq.gz | cut -c -11,14- | uniq)  ## Based on File name D3_S1_L001_R1_001.fastq.gz
do
var1=$(echo ${i}| cut -c -10) #Splitting the script to get the first part D3_S1_L001
var2=$(echo ${i}| cut -c 13-) #Splitting the string to get the second part _001.fastq.gz
##Checking the quality using fastqc. Note we have to change the directory of the programs as required
./fastqc $var1'_R1_'$var2 -o $var1'_R1_QC_'$var2
./fastqc $var1'_R2_'$var2 -o $var1'_R2_QC_'$var2

##Trimming required if we have the Phred score less than 30. Use trimming fuction as required/


done

##Printing  the end of the program
echo "end"

