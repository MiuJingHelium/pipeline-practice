#!/bin/bash

export LSF_DOCKER_VOLUMES="/storage1/fs1/martyomov/Active:/storage1/fs1/martyomov/Active /scratch1/fs1/martyomov:/scratch1/fs1/martyomov"
cd /scratch1/fs1/martyomov/carisa/GzmK_practice/test_argparseR


LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q general -G compute-martyomov \
     -J analysis.stg -n 5 -M 60000000 -o analysis.string.out \
     -e analysis.string.err -R 'select[mem>60000] rusage[mem=60000]' \
     -a "docker(psandhey/seurat_docker)" /bin/bash -c "Rscript analysis.R merge /scratch1/fs1/martyomov/carisa/GzmK_practice/results/align"
