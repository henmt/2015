#!/bin/bash

. sacgf_init.sh

load_module bowtie
load_module tophat

FASTQ1=${1}
FASTQ2=${2}
MATE_INNER_DIST=${3}
MATE_STD_DEV=${4:-20}

ORGANISM=${ORGANISM:-${DEFAULT_ORGANISM}}
BUILD=${BUILD:-${DEFAULT_BUILD}}

REFERENCE_GTF=$(get_igenomes_reference_gtf)
BOWTIE2_INDEX=$(get_igenomes_tool_index Bowtie2)
BOWTIE2_INDEX="${BOWTIE2_INDEX%.*}" # bowtie doesn't use .fa on the end

echo "Using reference genome ${REFERENCE_GENOME}"
echo "Using reference gtf ${REFERENCE_GTF}"

if [ $# -lt 3 ]; then
	echo "Usage $(basename $0) : read1.fastq read2.fastq innerdist (mate_std_dev)" >&2
	exit 1;
fi

if [ -z ${CORES} ]; then
	echo "No cores supplied - defaulting to using 4 cores."
	CORES=4
else
	echo "Using ${CORES} cores"
fi	

# Removed reference to inner dist etc in filename, as we can retrieve this from bam header
NAME=$(get_name ${FASTQ1})
NAME=${NAME%_R1}
TOPHAT_OUT_DIR=tophat2_out.${NAME}.innerdist_${MATE_INNER_DIST}.std_dev_${MATE_STD_DEV}.${BUILD}

tophat2 -p ${CORES} -o ${TOPHAT_OUT_DIR} --GTF ${REFERENCE_GTF} --mate-inner-dist=${MATE_INNER_DIST} --mate-std-dev=${MATE_STD_DEV} ${BOWTIE2_INDEX} ${FASTQ1} ${FASTQ2}

ALIGNER_ID=tophat2_pe

# Move & nicely rename files
for BED_TYPE in deletions insertions junctions; do
	mv ${TOPHAT_OUT_DIR}/${BED_TYPE}.bed ${NAME}.${BED_TYPE}.${ALIGNER_ID}.${BUILD}.bed
done

BAM_NAME=${NAME}.${ALIGNER_ID}.${BUILD}.bam
mv ${TOPHAT_OUT_DIR}/accepted_hits.bam ${BAM_NAME}
bam_index.sh ${BAM_NAME}

