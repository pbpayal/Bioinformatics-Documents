cd /Volumes/Drive\ 1_2/Payal/Genomes/Human/G38/Human_HISAT_Reference
for file in $(ls *sam)
do
  echo "My file name is $file"
  echo "Running ----SAMTOOLS SORT AND CONVERSION TO BAM----"
  samtools sort  -o sorted$file -O bam $file
  echo "The  output is sorted$file"
  
done

echo "end"