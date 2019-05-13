cut -f 1 counts_single_filtered.txt > var1
awk '$2 > 0' counts_single_filtered.txt > var2
for i in var1
	grep "ENSG00000280113" counts_paired.txt
	
