#!/bin/bash

cd /data/alignment_rescued2/

for f in $(ls *vdjca)
do
        mixcr assemble $f clones_$f.clns
done

exit
