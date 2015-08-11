#!/bin/bash

. sacgf_init.sh
load_module fastQC

FASTQ=${1}

# Fail on errors
set -e

if [ -z ${FASTQ} ]; then
        echo "No fastq file provided as argument"
        exit;
elif [ ! -e ${FASTQ} ]; then
        echo "No fastq file ('${FASTQ}') found"
        exit 1
else
        echo `basename ${0}` ${FASTQ}
fi

fastqc -o . --noextract ${FASTQ}
