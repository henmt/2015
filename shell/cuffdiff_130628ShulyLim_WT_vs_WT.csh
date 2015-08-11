#!/bin/csh

### Job name
#PBS -N cuffdiff_130628ShulyLim_f2_wt_vs_f5_wt

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
setenv CUFFDIFF_RUN_NAME f2_wt_vs_f5_wt
setenv WORKING_DIR /scratch/dlawrence/cuffdiff_${EXPERIMENT_NAME}

setenv ORGANISM Mus_musculus
setenv BUILD mm10

setenv CORES 6
setenv REFERENCE_GTF /data/sacgf/reference/iGenomes/${ORGANISM}/UCSC/${BUILD}/Annotation/Genes/genes.gtf

setenv BAMS_DIR /data/sacgf/molpath/data/aligned/130628ShulyLim
setenv F2_WILDTYPE ${BAMS_DIR}/WT-F2-Rep1_ATCACG_L007_R1.trimmed.tophat2_pe.mm10.bam,${BAMS_DIR}/WT-F2-Rep2_TAGCTT_L007_R1.trimmed.tophat2_pe.mm10.bam
setenv F5_WILDTYPE ${BAMS_DIR}/WT-F5-Rep1_ACTTGA_L007_R1.trimmed.tophat2_pe.mm10.bam,${BAMS_DIR}/WT-F5-Rep2_GGCTAC_L007_R1.trimmed.tophat2_pe.mm10.bam,${BAMS_DIR}/WT-F5-Rep3_GTGGCC_L007_R1.trimmed.tophat2_pe.mm10.bam

mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

cuffdiff -p ${CORES} -v ${REFERENCE_GTF} -L F2_WT,F5_WT -o ${CUFFDIFF_RUN_NAME} ${F2_WILDTYPE} ${F5_WILDTYPE}

