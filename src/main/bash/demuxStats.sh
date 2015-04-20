#!/bin/bash

mainDir=$1

for plate in `echo $mainDir*R1`
do
  out="$mainDir/$plate/counts.csv"
  all=`cat $mainDir/$plate/*/*R1*fastq|wc -l`
  invalid=`cat $mainDir/$plate/invalid/*R1*fastq|wc -l`
  undetermined=`cat $mainDir/$plate/undetermined/*R1*fastq|wc -l`

  invalidRate=`echo "scale=4;$invalid/$all" | bc`
  undeterminedRate=`echo "scale=4;$undetermined/$all" | bc`

  echo "invalid,undetermined">$out
  echo "$invalidRate,$undeterminedRate">>$out
done

out="$mainDir/read_counts.csv"
echo "plate,well,bc,nreads,nmapped">$out
for plate in `echo $mainDir*R1`
do
  for well in `ls $mainDir/$plate|grep -P "^[A-D]_|[A-H]\d{2}$"`
  do
    # echo "$plate/$well"
    if [[-f "$mainDir/$plate/$well/OR.fastq"]]
    then
      nreads=$((`wc -l $mainDir/$plate/$well/OR.fastq|cut -f1 -d" "` / 4))
      nmapped=`cut -f5 $mainDir/$plate/$well/OR.sam|tail -n+4|grep -v 0|wc -l`
      echo "$mainDir/$plate,$well,NA,$nreads,$nmapped">>$out
    else
      for bc in `ls OR_BC*.fastq|cat|cut -c6`
      do
        nreads=$((`wc -l $mainDir/$plate/$well/OR_BC$bc.fastq|cut -f1 -d" "` / 4))
        nmapped=`cut -f5 $mainDir/$plate/$well/OR_BC$bc.sam|tail -n+4|grep -v 0|wc -l`
        echo "$plate,$well,$bc,$nreads,$nmapped">>$out
      done
    fi
  done
done
