#!/bin/bash

find . -name "*.py" -not -path "*/1-advanced/*" -not -path "*/2-benchmark/*" | while read i
do
  echo $i
  grep -q "This example runs slow" $i
  if [ $? == 0 ]; then
    # If the example runs too long, skip the example
    echo "Skip slow example $i"
  else
    OMP_NUM_THREADS=2 python $i
  fi
done
