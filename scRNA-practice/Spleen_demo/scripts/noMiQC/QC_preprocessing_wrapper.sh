#!/bin/bash
indir=$1
outdir=$2
chmod +x ./QC_preprocessing_Single.sh
for sample in $(ls $indir); do
	/scratch1/fs1/martyomov/carisa/GzmK_practice/Spleen_demo/scripts/noMiQC/QC_preprocessing_Single.sh $indir $outdir $sample
done 

