#!/bin/bash

mainDir=$1

for plate in `echo $mainDir*R1`
do
	out="$plate/counts.csv"
	all=`cat $plate/*/*R1*fastq|wc -l`
	invalid=`cat $plate/invalid/*R1*fastq|wc -l`
	undetermined=`cat $plate/undetermined/*R1*fastq|wc -l`

	invalidRate=`echo "scale=4;$invalid/$all" | bc`
	undeterminedRate=`echo "scale=4;$undetermined/$all" | bc`

	echo "invalid,undetermined">$out
	echo "$invalidRate,$undeterminedRate">>$out
done
