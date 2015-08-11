INTRO
-----

To perform the analysis, we need to collect data from each bam, then combine this data together.

The bam processing can be done in parallel, so has been separated into its own script.

DESCRIPTION
-----------

-Measure read length distributions (all reads and piRNAs)

-Count primary and secondary piRNAs

-Count uridylation and adenylation


RUNNING
-------

* Step 0. Trim Fastqs, run FastQC on the trimmed file. Align the bams

* Step 1. Run the bam processing scripts.

# Your genome path will vary
GENOME_FASTA=/data/sacgf/reference/iGenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa

python lim_pirna_bam_processor.py --genome_fasta=${GENOME_FASTA} 

* Step 3. Run analysis script(s)

python lim_pirna_sequence_length_analysis.py ${FASTQC_DIR}

# This looks for the output from step 1 in the current directory. If this is from somewhere else, run:

python lim_pirna_analysis.py ${INPUT_DIR}

INSTALLATION NOTES
------------------

Developed on Ubuntu v13.04, Python v2.7.4

You may need to install some packages first. On my machine (Ubuntu v13.04) I have to run:

	sudo apt-get install g++ libpng-dev libfreetype6-dev python-dev python-setuptools python-pip gfortran libopenblas-dev liblapack-dev zlib1g-dev

Then you should be able to install the required packages by running:

	sudo pip install -U -r requirements.txt


