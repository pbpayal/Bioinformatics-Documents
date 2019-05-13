#!/bin/bash

cd /data/

for f in *R1.fastq
do
        r1=$f;
        r2=${f/R1.fastq/}R2.fastq;
        echo "The R1 file is $r1 and the R2 file is  $r2";
        mixcr align -p rna-seq -OallowPartialAlignments=true $r1 $r2 alignments_rna_$r1.vdjca
        echo "The output file is alignments_rna_$r1.vdjca"       
done

echo "Finished running mixcr"