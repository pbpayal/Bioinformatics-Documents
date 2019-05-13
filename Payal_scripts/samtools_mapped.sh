#!/bin/bash

cd /Volumes/Drive_1_2/Payal/Genmed/Val_Stephanie/Val\ Run\ at\ GU/CNMC-96662566 

for file in $(ls *bam)
do
samtools view -b -F 4 $file > mapped_$file
done

exit
