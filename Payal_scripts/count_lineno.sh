## Count the number of lines in a file without the first row or header row

cd /Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Cancer/96_Mike_CN_2021/convert_output 
for f in $(ls convert_output*)
do 
	awk 'END{print NR-1}' $f > clonotype_counts_$f
done
exit

