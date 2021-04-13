cd /data3/Payal/G176_Maria

#----------------------------------------TOPHAT---------------------------------------------
for file in $(ls *R1_001.fastq.gz)
do
  echo "My file name is $file"
  echo "Running ----Cutadapt-----"
  cutadapt -u 5 -o trimmed_$file $file
  echo "The output is trimmed_$file"
done
echo "------------------------------Finished Running Cutadapt------------------------------"

#----------------------------------------HISAT2---------------------------------------------
for file in $(ls *001.fastq.gz)
do
  echo "My file name is $file"
  echo "Running ----HISAT2-----"
  hisat2 --threads 8 --time -x /data3/Payal/Genomes/Mouse/GRCm38/Hisat2_Index/Mus_musculus.GRCm38.dna.primary_assembly -U $file -S aligned_$file.sam --summary-file summary_$file.txt
  echo "The output is aligned_$file.sam"
done
echo "------------------------------Finished Running HISAT2------------------------------"

#----------------------------------------Samtools Sort and Convert BAM---------------------------------------------
for file in $(ls aligned_*)
do
  echo "Running ----SAMTOOLS SORT AND CONVERSION TO BAM----"
  echo "My file name is $file"
  samtools sort -o sorted_$file -O bam $file
  echo "The output is sorted_$file"
  
done
echo "------------------------------Finished sorting and converting files with SAMTOOLS------------------------------"

#----------------------------------------Htseq Counts---------------------------------------------
for file in $(ls sorted_*) 
do
  echo "Running ----HTSEQ COUNTS----"
  echo "My file name is $file"
  htseq-count --format bam --order pos -t exon $file /data3/Payal/Genomes/Mouse/GRCm38/Mus_musculus.GRCm38.91.gtf > counts_$file.txt
  echo "The output is counts_$file.txt"
done
echo "------------------------------Finished HTSEQ Counts------------------------------"
echo "------------------------------Pipeline Over------------------------------"



tophat -p 8 -r 50 --bowtie1 -o tophat1  --no-novel-juncs --no-coverage-search -G /Volumes/Drive_1_2/Payal/Genomes/Human/G38/Homo_sapiens.GRCh38.91_ERCC92.gtf /Volumes/Drive_1_2/Payal/Genomes/Human/G38/Human_ERCC/Homo_sapiens.GRCh38.dna.toplevel_ERCC92 trimmed1_TS9_S27_L001_R1_001.fastq.gz TS9_S27_L001_R3_001.fastq.gz

tophat -p 8 -r 50 --bowtie1 --no-novel-juncs --no-coverage-search -G /Volumes/Drive_1_2/Payal/Genomes/Human/G38/Homo_sapiens.GRCh38.91_ERCC92.gtf /Volumes/Drive_1_2/Payal/Genomes/Human/G38/Human_ERCC/Homo_sapiens.GRCh38.dna.toplevel_ERCC92 TS9_S27_L001_R1_001.fastq.gz TS9_S27_L001_R3_001.fastq.gz
