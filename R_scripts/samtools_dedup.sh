#cd /Volumes/Drive 1_2/Payal/Kazue_Toru/G151
for file in $(ls sortedtrimmed*)
do
  echo "My file name is $file"
  echo "Running ----SAMTOOLS DEDUP----"
  samtools rmdup $file dedup$file
  echo "The  output is dedup$file"
  
done

echo "end"