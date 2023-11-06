#!/bin/bash

indir=$1

export LSF_DOCKER_VOLUMES="/storage1/fs1/martyomov/Active/:/storage1/fs1/martyomov/Active/ /scratch1/fs1/martyomov/carisa:/scratch1/fs1/martyomov/carisa /home/carisa:/home/carisa"

WD="$indir"
cd $WD
#current_time=$(date "+%Y.%m.%d-%H.%M.%S")
#export P_LOG="/scratch1/fs1/martyomov/carisa/logs/$current_time.log"

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q martyomov -G compute-martyomov \
        -n 2 -o Object_Creation_NKILC.out \
        -e Object_Creation_NKILC.err -R 'rusage[mem=32GB] span[hosts=1]' \
	-J Object_Creation_NKILC -M 32GB \
        -a "docker(kalisaz/scrna-extra:r4.3.0)" /bin/bash -c \
	"Rscript ./Object_Creation_NKILC.R"
