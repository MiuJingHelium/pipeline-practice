#!/bin/bash

for dir in */; do
	files=$(ls $dir)
	sample=${dir%/*}
	echo "$sample"
	for f in $files; do
		type=${f#*$sample}
		type=${type#*_}
		#echo "$type"
		mv $dir$f $dir$type
	done
done
