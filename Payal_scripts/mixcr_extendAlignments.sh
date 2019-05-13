#!/bin/bash

cd /data/alignment_rescued2/

for f in $(ls *vdjca)
do
	mixcr extendAlignments $f extended_$f 
done

exit

