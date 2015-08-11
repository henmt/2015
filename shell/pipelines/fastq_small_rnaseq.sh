#!/bin/bash
# same as fastq_small_mirna_count but without the count (ie not keeping SAM)

. sacgf_init.sh

FASTQ=$1
BUILD=${BUILD:-${DEFAULT_BUILD}}

if [ -z ${FASTQ} ]; then
        echo "No fastq file provided as argument"
        exit 1;
elif [ ! -e ${FASTQ} ]; then
        echo "No fastq file ('${FASTQ}') found"
        exit 1;
else
        echo "FASTQ = '${FASTQ}'";
fi

NAME=$(get_name ${FASTQ})
EXTENSION=$(get_extension ${FASTQ})

#check whether output file exists
NOLINKERS=${NAME}.cutadapt_trimmed
echo "Checking for '${NOLINKERS}.${EXTENSION}'" 1>&2
if [ -e ${NOLINKERS}.${EXTENSION} ]; then
	echo "${NOLINKERS}.${EXTENSION} already exists." 1>&2
else
	# Will guess adapter
	CUTADAPT_OUTPUT=${NAME}.cutadapt.output.txt
	echo "Storing cutadapt output to '${CUTADAPT_OUTPUT}'"
	fastx_trim_adapter.sh ${FASTQ} > ${CUTADAPT_OUTPUT}
fi

#make a fastqc report if one doesnt exist already
FASTQC_ZIP=${NOLINKERS}_fastqc.zip
if [ ! -e ${FASTQC_ZIP} ]; then
        fastqc.sh ${NOLINKERS}.${EXTENSION}
fi

#align using bwa
ALIGNED_NAME=${NOLINKERS}.bwa.${BUILD}
BAM=${ALIGNED_NAME}.bam
if [ ! -e ${BAM} ]; then
	echo "BWA aligning ${SAM}"
	fastq_bwa_to_bam.sh ${NOLINKERS}.${EXTENSION}
fi
