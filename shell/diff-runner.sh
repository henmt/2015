#!/bin/sh

# Unfiltered changes.
GOMAXPROCS=4 parallel -j 2 -N 3 ~/Development/src/bitbucket.org/sacgf/lim_2014/go/length-heat-annot-diff -ref mm10.mfa -in {2},{3} -f {1} -out filtered/{3.}-diff ::: 0 {wt,mut}-f2.bam 1 {wt,mut}-f2.bam 2 {wt,mut}-f2.bam 0 {wt,mut}-f5.bam 1 {wt,mut}-f5.bam 2 {wt,mut}-f5.bam

# Changes filtered for a variety of features.
for CLASS in repeat/{SINE,LINE,LTR} CDS exon intron
do
	GOMAXPROCS=4 parallel -j 2 -N 3 ~/Development/src/bitbucket.org/sacgf/lim_2014/go/length-heat-annot-diff -ref mm10.mfa -in {2},{3} -f {1} -annot mm10.feat.gff -class $CLASS -out filtered/{3.}-${CLASS##repeat/}-diff ::: 0 {wt,mut}-f2.bam 1 {wt,mut}-f2.bam 2 {wt,mut}-f2.bam 0 {wt,mut}-f5.bam 1 {wt,mut}-f5.bam 2 {wt,mut}-f5.bam
done
