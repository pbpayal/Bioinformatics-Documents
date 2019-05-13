#!/bin/bash

cd /Users/pbanerjee/Desktop/Katrina_counts/Cell_Line1

for file in $(ls counts_star*)
do 
	head -n $(( $(wc -l $file | awk '{print $1}') - 5 )) $file > extr_$file
done

# But this creates a blank line at the end of the file
	
