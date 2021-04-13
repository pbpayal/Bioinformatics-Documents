
cd 

for file in *R1_001.fastq.gz
do 
mv $file ${file//out_fwd_paired_/}
done