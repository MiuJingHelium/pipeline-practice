#!/bin/bash

indir=$1

export LSF_DOCKER_VOLUMES="/storage1/fs1/martyomov/Active/:/storage1/fs1/martyomov/Active/ /scratch1/fs1/martyomov/carisa:/scratch1/fs1/martyomov/carisa /home/carisa:/home/carisa"

WD="$indir"
cd $WD

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q martyomov -G compute-martyomov \
        -n 2 -o Check_Object_CD4_TS.out \
        -e Check_Object_CD4_TS.err -R 'rusage[mem=32GB] span[hosts=1]' \
	-J Check_Object_CD4_TS -M 32GB \
        -a "docker(kalisaz/scrna-extra:r4.3.0)" /bin/bash -c \
	"Rscript ./Check_Object_CD4_TS.R"
