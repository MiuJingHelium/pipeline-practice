#!/bin/bash

INDIR=$1
OUTDIR=$2
REF='/storage1/fs1/martyomov/Active/References/10X/SC/Mouse/refdata-cellranger-mm10-3.0.0'
#WD=`pwd`
for SAMPLE in $(ls $INDIR);do
	#echo "align_task.bash $INDIR $OUTDIR $SAMPLE $REF"
	#mkdir -p $OUTDIR/$SAMPLE
	/scratch1/fs1/martyomov/carisa/GzmK_practice/Spleen_demo/scripts/align_task.bash $INDIR $OUTDIR $SAMPLE $REF
done
#$4 transcriptome
