cd /Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Feature_counts2
for file in $(ls *txt)
do 
	head -n 53465 $file | awk '{print $1,$7}' | awk 'NR>2' > deseq2_$file
done

exit
