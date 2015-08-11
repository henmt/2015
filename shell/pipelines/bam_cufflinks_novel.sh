#!/bin/bash

. sacgf_init.sh
load_module cufflinks

# Setup
if [ $# -ne 1 ]; then
	echo "Usage $(basename $0): tophat.bam"
	exit 1
fi
BAM=$1

ORGANISM=${ORGANISM:-${DEFAULT_ORGANISM}}
BUILD=${BUILD:-${DEFAULT_BUILD}}
REFERENCE_GTF=$(get_igenomes_reference_gtf)

if [ -z ${CORES} ]; then
	echo "No cores supplied - defaulting to using 4 cores."
	CORES=4
else
	echo "Using ${CORES} cores"
fi	

BAMNAME=$(basename ${BAM} .bam)
BAM_BUILD=$(get_extension ${BAMNAME})
if [ ${BUILD} != ${BAM_BUILD} ]; then
	echo "BUILD (${BUILD}) != BAM BUILD (${BAM_BUILD})" >&2
	exit 1
fi

CUFFLINKS_DIR=cufflinks_${BAMNAME}.${BUILD}
cufflinks -p ${CORES} --quiet -g ${REFERENCE_GTF} -o ${CUFFLINKS_DIR} ${BAM} #-g = also assemble novel

