cd /Volumes/Drive\ 1_2/Payal/Kazue_Toru/G151
for file in $(ls sortedtrimmed*)
do
  echo "My file name is $file"
  echo "Removing ----CHR FROM BAM----"
  samtools view -h $file | sed 's/chr//' | samtools view -Shb - -o chr_rem$file
  echo "The  output is chr_rem$file"
  
done

echo "end"
