#!/bin/csh

#PBS -N SL_WT_F5_Rep1_R1.fastq

#PBS -j oe

#PBS -m ae

#PBS -M david.lawrence@health.sa.gov.au


#PBS -l ncpus=4
#PBS -l mem=12gb,vmem=12gb
#PBS -l walltime=400:00:00

setenv ORGANISM Mus_musculus
setenv BUILD mm10
setenv WORKDIR /scratch/dlawrence/fastq_small_rnaseq.sh
setenv CORES 4
setenv PATH /home/users/dlawrence/localwork/bioinformatics/scripts/intercept:${PATH}:/home/users/dlawrence/localwork/bioinformatics/scripts/bin:/home/users/dlawrence/localwork/bioinformatics/scripts/pipelines:/home/users/dlawrence/localwork/bioinformatics/scripts/tasks

mkdir -p ${WORKDIR}
cd ${WORKDIR}

fastq_small_rnaseq.sh /data/sacgf/molpath/data/unaligned/130517ShulyLim/fastq/SL_WT_F5_Rep1_R1.fastq.gz
