#cd /Volumes/Drive 1_2/Payal/Kazue_Toru/G151
for file in $(ls *.fastq.gz.bam)
do
  echo "My file name is $file"
  echo "Running ----FLAGSTAT----"
  samtools flagstat --threads 2 $file > flagstat$file
  echo "The  output is flagstat$file"
  
done

echo "end"