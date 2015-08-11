#!/bin/bash

# Quality Trim
# Tophat
# Cuffdiff w/Novel

. sacgf_init.sh
load_module cutadapt

FASTQ1=${1}
FASTQ2=${2}

ORGANISM=${ORGANISM:-${DEFAULT_ORGANISM}}
BUILD=${BUILD:-${DEFAULT_BUILD}}

MIN_LENGTH=18
QUALITY_CUTOFF=28

echo "Using build ${BUILD}"

if [ -z ${FASTQ1} ] || [ -z ${FASTQ2} ]; then
	echo "Usage $(basename $0): read1.fastq read2.fastq (innerdist) (mate_std_dev)"
	exit 1;
fi

function trim_and_fastqc {
	FASTQ=$1
	TRIMMED_FASTQ=$2

	cutadapt --quality-cutoff=${QUALITY_CUTOFF} --minimum-length=${MIN_LENGTH} --output=${TRIMMED_FASTQ} ${FASTQ} > ${TRIMMED_FASTQ}.cutadapt.log
	fastqc.sh ${TRIMMED_FASTQ}
}


echo "Trimming 2 fastq files in separate threads"
TRIMMED_FASTQ1=$(get_name ${FASTQ1}).trimmed.$(get_extension ${FASTQ1})
trim_and_fastqc ${FASTQ1} ${TRIMMED_FASTQ1} &

TRIMMED_FASTQ2=$(get_name ${FASTQ2}).trimmed.$(get_extension ${FASTQ2})
trim_and_fastqc ${FASTQ2} ${TRIMMED_FASTQ2} &
wait

echo "Both files trimmed"

NAME=$(get_name ${TRIMMED_FASTQ1})
NAME=${NAME%_R1}
ALIGNER_ID=tophat2_pe
TOPHAT_BAM=${NAME}.${ALIGNER_ID}.${BUILD}.bam

echo "Tophat"
fastq_tophat_pe.sh ${TRIMMED_FASTQ1} ${TRIMMED_FASTQ2} ${3} ${4}

echo "Cufflinks Novel"
bam_cufflinks_novel.sh ${TOPHAT_BAM}

