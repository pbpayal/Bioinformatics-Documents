GATK Test

** Reference Genome - hg38 from ENSEMBL - FASTA file
* Create index file - samtool faidx
* Create .dict file - picard CreateSequenceDictionary R=input.fa O=output.dict


* Download NA12878 CRAM file from 
* Convert CRAM file to BAM file - Samtools 
* Markduplicates - Picard with Remove duplicates

