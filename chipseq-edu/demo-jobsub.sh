#!/bin/bash

# Launch in compute client node (not in interactive job session):

export STORAGE1="/storage1/fs1/martyomov/Active"
export SCRATCH1="/scratch1/fs1/martyomov"
export LSF_DOCKER_VOLUMES="$STORAGE1:$STORAGE1 $SCRATCH1:$SCRATCH1 $HOME:$HOME"

cd /scratch1/fs1/martyomov/carisa/chipseq-smk-pipeline-edu
P_WD=`pwd`
LSF_DOCKER_ENV_FILE=$P_WD/lsf_docker_env_file.env
export SMK_DOCKER_IMG="snakemake/snakemake:stable"
export P_LOG=$P_WD/logs/pipeline.log
export PATH="/opt/conda/envs/snakemake/bin:/opt/conda/bin:$PATH" 

# * optionally you could add `--restart-times 0` to override `lsf_demo` 
#    profile setting
L_CORES=4; LSF_DOCKER_ENV_FILE=$P_WD/lsf_docker_env_file.env; mkdir -p logs; bsub -cwd $HOME -n $L_CORES -G compute-martyomov -q general -oo $P_LOG -R 'span[hosts=1]' -a "docker($SMK_DOCKER_IMG)" /usr/bin/script -fqe /dev/null  -c "source /etc/bash.bashrc; cd $P_WD; export TMPDIR=$P_WD/tmp; snakemake -j 50 --profile lsf_demo --local-cores $L_CORES"

# track snakemake logs
touch -f $P_LOG && tail -f $P_LOG
 
# for colored logs don't forget `-R` in less:
less -R $P_LOG
