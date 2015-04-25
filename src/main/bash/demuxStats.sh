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

out="${mainDir}read_counts.csv"
echo "plate,well,bc,nreads,nmapped">$out
for plate in `echo $mainDir*R1`
do
  for well in `ls $plate|grep -P "^[A-D]_|[A-H]\d{2}$"`
  do
    # echo "$plate/$well"
    if [[ -f "$plate/$well/OR.fastq" ]]
    then
      nreads=$((`wc -l $plate/$well/OR.fastq|cut -f1 -d" "` / 4))
      nmapped=`cut -f5 $plate/$well/OR.sam|tail -n+4|grep -v 0|wc -l`
      echo "$plate,$well,NA,$nreads,$nmapped">>$out
    else
      for bc in `ls OR_BC*.fastq|cat|cut -c6`
      do
        nreads=$((`wc -l $plate/$well/OR_BC$bc.fastq|cut -f1 -d" "` / 4))
        nmapped=`cut -f5 $plate/$well/OR_BC$bc.sam|tail -n+4|grep -v 0|wc -l`
        echo "$plate,$well,$bc,$nreads,$nmapped">>$out
      done
    fi
  done
done
