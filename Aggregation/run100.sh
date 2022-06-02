#!/bin/bash
SEED=400;

for((i=101;i<=$SEED;i++));
do
    seed=$i
    echo $seed #this prints the seed
    nohup ./main100.out $seed > ./outputs/myscreenout100_$seed
done
