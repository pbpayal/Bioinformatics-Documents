#!/bin/bash

cd /data3/Payal/Projects/G216_Toru/Fastq_Files

for f in *1_001_trimmed.fastq
do
        r1=$f;
        r2=${f/1_001_trimmed.fastq/}3_001.fastq;
        echo "The R1 file is $r1 and the R2 file is  $r2";
        hisat2 --threads 8 -t --summary-file summary_$f.txt -p 8 -x /data3/Payal/Genomes/ERCC_Genome/Hisat2_ERCC_Index/ERCC92  -1 $r1 -2 $r2 -S aligned_$r1.sam
        echo "The output file is aligned_$r1.sam"
done

echo "Finished running Hisat2"
