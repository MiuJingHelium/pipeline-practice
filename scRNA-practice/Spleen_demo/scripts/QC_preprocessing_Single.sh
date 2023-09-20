#!/bin/bash

sample=$3
outdir=$2
indir=$1
export LSF_DOCKER_VOLUMES="/storage1/fs1/martyomov/Active/:/storage1/fs1/martyomov/Active/ /scratch1/fs1/martyomov/carisa:/scratch1/fs1/martyomov/carisa /home/carisa:/home/carisa"
mkdir -p "$outdir/$sample"
WD="$outdir$sample"
cd $WD
#current_time=$(date "+%Y.%m.%d-%H.%M.%S")
#export P_LOG="/scratch1/fs1/martyomov/carisa/logs/$current_time.log"

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q martyomov -G compute-martyomov \
        -n 2 -o QCPrep_single.out \
        -e QCPrep_single.err -R 'rusage[mem=10GB] span[hosts=1]' \
	-J QCPrep_single -M 10GB \
        -a "docker(kalisaz/scrna_prep:test)" /bin/bash -c \
	"Rscript /scratch1/fs1/martyomov/carisa/GzmK_practice/Spleen_demo/scripts/QC_preprocessing_Single.R $indir $outdir$sample $sample"
