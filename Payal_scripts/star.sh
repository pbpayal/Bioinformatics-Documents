#!/bin/bash

cd /data3/Payal/Projects/Evan_G276/trimmed_paired_fastq_files


for f in *R1.fastq.gz
do
	r1=$f;
	r2=${f/R1.fastq.gz/}R2.fastq.gz  ;
	echo "The R1 file is $r1 and the R2 file is $r2";
	STAR --runThreadN 8 --runMode alignReads --limitGenomeGenerateRAM=119000000000 --genomeSAsparseD 3 --genomeSAindexNbases 12 --genomeChrBinNbits=16 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix star_out_$r1 --genomeDir /data3/Payal/Genomes/Mouse/GRCm38_EGFP/STAR_GRCm38_EGFP/  --sjdbGTFfile /data3/Payal/Genomes/Mouse/GRCm38_EGFP/Mus_musculus.GRCm38.91_EGFP.gtf --readFilesIn $r1 $r2
	echo "The output file is star_out_$r1"	
done

echo "Finished running STAR"