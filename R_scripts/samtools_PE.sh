#!/bin/bash

cd /Volumes/Drive_1_2/Payal/Genmed/Val_Stephanie/Val\ Run\ at\ GU/CNMC-96662566

for file in $(ls mapped*)
do
samtools view -bf 1 $file > PE_$file.bam
done

exit
