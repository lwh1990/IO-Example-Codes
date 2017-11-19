#!/bin/bash
FILES=`find . | grep job.sh`
ACCOUNT=$1
for f in $FILES
do
	sed -i "s/#PBS -A.*$/#PBS -A ${ACCOUNT}/" $f
	echo "Account ${ACCOUNT} set for ${f}"
done
