#!/bin/bash

for f in $(ls *vdjca)
do
	mixcr assemblePartial $f alignmentsRescued1_$f
done

exit
