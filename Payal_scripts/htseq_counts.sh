cd /Volumes/Drive\ 1_2/Payal/Kazue_Toru/G151
for file in $(ls chr_remsortedtrimmed*)
do
  echo "My file name is $file"
  echo "Running ----HTSEQ COUNTS----"
  htseq-count --format bam --order pos -t exon $file Mus_musculus.GRCm38.91.gtf > counts$file.txt
  echo "The  output is counts$file.txt"
  
done

echo "end"