#!/bin/csh

### Job name
#PBS -N cuffdiff_130628ShulyLim_F2_reference_only

### Join queuing system output and error files into a single output file
#PBS -j oe

### Send email to user when job ends or aborts
#PBS -m ae

### email address for user
#PBS -M david.lawrence@health.sa.gov.au

### Queue name that job is submitted to
###PBS -q sacgf

### Request nodes, memory, walltime. NB THESE ARE REQUIRED
#PBS -l ncpus=6
#PBS -l mem=15gb,vmem=15gb
#PBS -l walltime=50:00:00

module load cufflinks

setenv EXPERIMENT_NAME 130628ShulyLim
setenv CUFFDIFF_RUN_NAME F2_reference_only
setenv WORKING_DIR /scratch/dlawrence/cuffdiff_${EXPERIMENT_NAME}

setenv ORGANISM Mus_musculus
setenv BUILD mm10

setenv CORES 6
setenv REFERENCE_GTF /data/sacgf/reference/iGenomes/${ORGANISM}/UCSC/${BUILD}/Annotation/Genes/genes.gtf

setenv BAMS_DIR /data/sacgf/molpath/data/aligned/130628ShulyLim
setenv WILDTYPE ${BAMS_DIR}/WT-F2-Rep1_ATCACG_L007_R1.trimmed.tophat2_pe.mm10.bam,${BAMS_DIR}/WT-F2-Rep2_TAGCTT_L007_R1.trimmed.tophat2_pe.mm10.bam
setenv MUTANT ${BAMS_DIR}/Mut-F2-Rep1_CGTACG_L007_R1.trimmed.tophat2_pe.mm10.bam,${BAMS_DIR}/Mut-F2-Rep2_GCCAAT_L008_R1.trimmed.tophat2_pe.mm10.bam,${BAMS_DIR}/Mut-F2-Rep3_GTGAAA_L008_R1.trimmed.tophat2_pe.mm10.bam

mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

cuffdiff -p ${CORES} -v ${REFERENCE_GTF} -L WT,Mutant -o ${CUFFDIFF_RUN_NAME} ${WILDTYPE} ${MUTANT}

