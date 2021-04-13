cd /Volumes/Drive\ 1_2/Payal/Kazue_Toru/G151

#----------------------------------------TOPHAT---------------------------------------------
for file in $(ls *_001.fastq.gz)
do
  echo "My file name is $file"
  echo "Running ----HISAT2-----"
  hisat2 -q --threads 4 --time -x genome -U $file -S aligned$file.sam
  echo "The  output is aligned$file.sam"
done

echo "------------------------------Finished Running TOPHAT------------------------------"

for file in $(ls *tophat_out)
	cd $file
	for file in $(ls *bam)
	do
  		echo "My file name is $file"
  		echo "Running ----SAMTOOLS SORT AND CONVERSION TO BAM----"
  		samtools sort  -o sorted$file -O bam $file
  		echo "The  output is sorted$file"
  
	done

	echo "------------------------------Finished sorting files with SAMTOOLS------------------------------"


	do
  		echo "My file name is $file"
  		echo "Running ----HTSEQ COUNTS----"
  		htseq-count --format bam --order pos -t exon $file Mus_musculus.GRCm38.91.gtf > counts$file.txt
  		echo "The  output is counts$file.txt"
  
	done

	echo "------------------------------Finished HTSEQ Counts------------------------------"