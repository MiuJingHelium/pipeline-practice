#!/bin/bash



INDIR=$1
OUTDIR=$2
NAME=$3
REF=$4

export LSF_DOCKER_VOLUMES="/storage1/fs1/martyomov/Active/:/storage1/fs1/martyomov/Active/  /scratch1/fs1/martyomov:/scratch1/fs1/martyomov /home/carisa:/home/carisa"
cd $OUTDIR


LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q general -G compute-martyomov \
        -J ${NAME}.stg -n 8 -M 64GB -o ${NAME}.align.out \
        -e ${NAME}.align.err -R 'select[mem>64MB] rusage[mem=64GB] span[hosts=1]' \
        -a "docker(psandhey/cellranger)" /bin/bash -c "cellranger count --id=$NAME \
        --transcriptome=/storage1/fs1/martyomov/Active/References/10X/SC/Mouse/refdata-cellranger-mm10-3.0.0 \
        --fastqs=${INDIR}/${NAME} \
        --localmem=64 \
        --localcores=16"
