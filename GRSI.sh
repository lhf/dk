#!/bin/sh

N=6

for R in `seq 0 $N`
do
	make dk R=$R O=dk-$R stats
done

for R in `seq 0 $N`
do
	make kf R=$R O=kf-$R
done

