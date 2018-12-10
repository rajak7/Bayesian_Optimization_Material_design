#!/bin/sh

for M1 in Mo W
do
  for X1 in S Se Te
  do
    for M2 in Mo W
    do
      for X2 in S Se Te
      do
        echo $M1${X1}2-$M2${X2}2
        python snl_prep.py -s $M1${X1}2-$M2${X2}2 -d .
      done
    done
  done
done
