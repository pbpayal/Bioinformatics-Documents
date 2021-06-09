cd /Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/James/James_htseq_reverse

for file in $(ls counts*)
do 
	head -n 53465 $file | awk '{print $1,$2}'  > deseq2_$file
done

exit