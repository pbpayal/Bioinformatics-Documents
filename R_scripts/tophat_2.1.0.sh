cd /Volumes/Drive\ 1_2/Payal/Kazue_Toru/G151

for file in $(ls *_001.fastq.gz)
do
  echo "My file name is $file"
  echo "Running ----TOPHAT 2.1.0-----"
  tophat -p 6  --bowtie1 -o tophat_out$file -r 100 --no-novel-juncs --no-coverage-search -G Mus_musculus.GRCm38.91.gtf Mus_musculus.GRCm38.dna.toplevel $file
  echo "The  output is trimmed$file.sam"
done

echo "end"
