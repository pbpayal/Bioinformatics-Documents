cd /Volumes/Drive\ 1_2/Payal/Kazue_Toru/G151
for file in $(ls *_001.fastq.gz)
do
  echo "My file name is $file"
  echo "Running ----HISAT2-----"
  hisat2 -q --threads 4 --time -x genome -U $file -S aligned$file.sam
  echo "The  output is trimmed$file.sam"
done

echo "end"
