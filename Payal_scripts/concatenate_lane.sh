#!/bin/bash

cd /Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Cancer/96_Mike_CN_2021/Old_Data_2018/Fastq_files

for f in $(ls *L001_R1_001.fastq.gz)
do
	r1=$f;
	r2=${f/L001_R1_001.fastq.gz/}L002_R1_001.fastq.gz;
	r3=${f/L001_R1_001.fastq.gz/}L003_R1_001.fastq.gz;
	r4=${f/L001_R1_001.fastq.gz/}L004_R1_001.fastq.gz;
	echo "The R1 file is $r1, R2 file is  $r2, R3 is $r3, R4 is $r4";
	cat $r1 $r2 $r3 $r4 > final_merged_$r1
	echo "The final filename is $final_merged_$r1"
done

for f in $(ls *L001_R2_001.fastq.gz)
do
	r1=$f;
	r2=${f/L001_R2_001.fastq.gz/}L002_R2_001.fastq.gz;
	r3=${f/L001_R2_001.fastq.gz/}L003_R2_001.fastq.gz;
	r4=${f/L001_R2_001.fastq.gz/}L004_R2_001.fastq.gz;
	echo "The R1 file is $r1, R2 file is  $r2, R3 is $r3, R4 is $r4";
	cat $r1 $r2 $r3 $r4 > final_merged_$r1
	echo "The final filename is $final_merged_$r1"
done
exit