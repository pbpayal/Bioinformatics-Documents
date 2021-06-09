#files -https://figshare.com/articles/Data_Carpentry_Genomics_beta_2_0/7726454


mkdir -p dc_workshop/raw_reads/../docs/../metadata/../qc/../trimmed_reads/../results/bam/../sam/../bcf/../vcf

cd dc_workshop

wget https://raw.githubusercontent.com/datacarpentry/wrangling-genomics/gh-pages/files/Ecoli_metadata_composite.csv

mv Ecoli_metadata_composite.csv metadata


cp ~/.miniconda3/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/NexteraPE-PE.fa metadata/

cp ~/.backup/untrimmed_fastq/*fastq.gz raw_reads/

fastqc -o qc/ raw_reads/*

for result in qc/*zip; do 
 unzip -d qc $result
 done
 

cat qc/*/summary.txt > docs/fastqc_summaries.txt


grep FAIL docs/fastqc_summaries.txt


for infile in raw_reads/*_1.fastq.gz
    do base=$(basename ${infile} _1.fastq.gz)
    trimmomatic PE ${infile} raw_reads/${base}_2.fastq.gz \
        trimmed_reads/${base}_1.trim.fastq.gz trimmed_reads/${base}_1un.trim.fastq.gz \
        trimmed_reads/${base}_2.trim.fastq.gz trimmed_reads/${base}_2un.trim.fastq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:metadata/NexteraPE-PE.fa:2:40:15 
done


#################################################################


mkdir dc_workshop/ref_genome

curl -L -o ref_genome/ecoli_rel606.fasta.gz \
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz

gunzip ref_genome/ecoli_rel606.fasta.gz


curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
tar xvf sub.tar.gz
mv sub/ trimmed_fastq_small


bwa index ref_genome/ecoli_rel606.fasta

bwa mem \
    ref_genome/ecoli_rel606.fasta \
    trimmed_fastq_small/SRR2584866_1.trim.sub.fastq \
    trimmed_fastq_small/SRR2584866_2.trim.sub.fastq \
    > results/sam/SRR2584866.aligned.sam

samtools view -S -b results/sam/SRR2584866.aligned.sam > results/bam/SRR2584866.aligned.bam

samtools index results/bam/SRR2584866.aligned.sorted.bam
samtools tview results/bam/SRR2584866.aligned.sorted.bam data/ref_genome/ecoli_rel606.fasta

samtools sort -o results/bam/SRR2584866.aligned.sorted.bam results/bam/SRR2584866.aligned.bam 

#bam file summary statistics
samtools flagstat results/bam/SRR2584866.aligned.sorted.bam

bcftools mpileup -O b -o results/bcf/SRR2584866_raw.bcf \
    -f ref_genome/ecoli_rel606.fasta \
    results/bam/SRR2584866.aligned.sorted.bam 

bcftools call --ploidy 1 -m -v \
    -o results/bcf/SRR2584866_variants.vcf \
    results/bcf/SRR2584866_raw.bcf 

bcftools call --ploidy 1 -m -v \
    -o results/bcf/SRR2584866_variants.vcf \
    results/bcf/SRR2584866_raw.bcf 

vcfutils.pl varFilter results/bcf/SRR2584866_variants.vcf > results/vcf/SRR2584866_final_variants.vcf






for fastq in raw_reads/*; do 
        echo $fastq
done

for fastq in raw_reads/SRR2584863_1.fastq.gz raw_reads/SRR2584863_2.fastq.gz raw_reads/SRR2584866_1.fastq.gz raw_reads/SRR2584866_2.fastq.gz raw_reads/SRR2589044_1.fastq.gz raw_reads/SRR2589044_2.fastq.gz; do
        echo $fastq; done


###TMUX

tmux new -s alignment
tmux attach -t alignment
tmux list-sessions



## CSV while
while read line; do echo "${line}; done < file.csv

while read line; do echo ${line} | do cut -f 8 -d "," ${line}; done < file.csv

while read line; do sample_id=$(echo ${line} | cut -f 8 -d ","); cit=$(echo ${line} |cut -f 12 -d ","); echo "${cit} ${sample_id}"; done < file.csv

