#!/bin/bash

cd /data3/Payal/G176_Maria/

for file in $(ls *.fastq.gz)
do
  echo "My file name is $file"
  echo "Running fastqc....."
  fastqc $file
done

echo "Saving FastQC Results....."
mv *fastqc.zip /data3/Payal/G176_Maria/FASTQC_Results_Zip/
mv *.html /data3/Payal/G176_Maria/FASTQC_Results_HTML/

cd /data3/Payal/G176_Maria/FASTQC_Results_Zip/

echo "Unzipping......"
for zip in *.zip
do 
	unzip $zip
done
