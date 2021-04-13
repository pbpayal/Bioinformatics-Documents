#!/bin/bash

cd /data/clones/

for f in $(ls *clns)
do
	mixcr  exportClones $f $f.txt
done

exit
