#!/bin/bash
indir=$1
outdir=$2
chmod +x ./QC_preprocessing_Single.sh
for sample in $(ls $indir); do
	./QC_preprocessing_Single.sh $indir $outdir $sample
done 

