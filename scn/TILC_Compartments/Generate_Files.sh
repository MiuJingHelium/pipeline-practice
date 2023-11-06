#!/bin/bash

indir=$1
#object=$2
#outdir=$2

export LSF_DOCKER_VOLUMES="/storage1/fs1/martyomov/Active/:/storage1/fs1/martyomov/Active/ /scratch1/fs1/martyomov/carisa:/scratch1/fs1/martyomov/carisa /home/carisa:/home/carisa"

WD="$indir"
cd $WD
#current_time=$(date "+%Y.%m.%d-%H.%M.%S")
#export P_LOG="/scratch1/fs1/martyomov/carisa/logs/$current_time.log"

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q martyomov -G compute-martyomov \
        -n 2 -o Generate_Files.out \
        -e Generate_Files.err -R 'rusage[mem=24GB] span[hosts=1]' \
	-J Generate_Files -M 24GB \
        -a "docker(kalisaz/scrna-extra:r4.3.0)" /bin/bash -c \
	"Rscript ./Generate_Files.R"
