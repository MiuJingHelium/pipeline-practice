#!/bin/bash

indir=$1 #An abosulte path

export LSF_DOCKER_VOLUMES="/storage1/fs1/martyomov/Active/:/storage1/fs1/martyomov/Active/ /scratch1/fs1/martyomov/carisa:/scratch1/fs1/martyomov/carisa /home/carisa:/home/carisa"

cd $indir
mkdir -p "$indir/DE"

#current_time=$(date "+%Y.%m.%d-%H.%M.%S")
#export P_LOG="/scratch1/fs1/martyomov/carisa/logs/$current_time.log"

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q martyomov -G compute-martyomov \
        -n 2 -o plot.out \
        -e plot.err -R 'rusage[mem=8GB] span[hosts=1]' \
	-J plot -M 8GB \
        -a "docker(kalisaz/scrna_analysis)" /bin/bash -c \
	"Rscript /scratch1/fs1/martyomov/carisa/GzmK_practice/Spleen_demo/scripts/Volcano_plot.R $indir"
