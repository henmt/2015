#!/bin/bash

# generate expression level data
GOMAXPROCS=4 parallel -j4 length-heat -in {} -f 1 -out heat-data/{.} ::: {wt,mut}-f2.bam {wt,mut}-f5.bam


# generate diff data for circle figures
GOMAXPROCS=4 parallel -j1 -N2 length-heat-annot-diff -denest -in {1},{2} -f 1 -out piRNA-expression/filtered-diffs/{2.}-diff ::: {wt,mut}-f2.bam {wt,mut}-f5.bam
for CLASS in repeat/{SINE,LINE,LTR} CDS exon intron
do
	GOMAXPROCS=4 parallel -j1 -N2 length-heat-annot-diff -denest -in {1},{2} -f 1 -annot mm10.feat.gff -class $CLASS -out piRNA-expression/filtered-diffs/{2.}-${CLASS##repeat/}-diff ::: {wt,mut}-f2.bam {wt,mut}-f5.bam
done

# generate diff count data for circle figures
GOMAXPROCS=4 parallel -j1 -N2 length-heat-annot-diff -denest -in {1},{2} -f 1 -out piRNA-expression/filtered-diff-counts/{2.}-diff ::: {wt,mut}-f2.bam {wt,mut}-f5.bam
for CLASS in repeat/{SINE,LINE,LTR} CDS exon intron
do
	GOMAXPROCS=4 parallel -j1 -N2 length-heat-annot-diff -denest -in {1},{2} -f 1 -annot mm10.feat.gff -class $CLASS -out piRNA-expression/filtered-diff-counts/{2.}-${CLASS##repeat/}-diff ::: {wt,mut}-f2.bam {wt,mut}-f5.bam
done

# generate data for end overlap analysis
GOMAXPROCS=6 parallel -j1 -N5 overlap-end -f {1} -pair {2} -pair {3} -pair {4} -out overlap-data/{5} ::: 0 SL_Mut_F2_Rep*bam mut-f2-reps 1 SL_Mut_F2_Rep*bam mut-f2-reps 2 SL_Mut_F2_Rep*bam mut-f2-reps 0 SL_Mut_F5_Rep*bam mut-f5-reps 1 SL_Mut_F5_Rep*bam mut-f5-reps 2 SL_Mut_F5_Rep*bam mut-f5-reps 0 SL_WT_F5_Rep*bam wt-f5-reps 1 SL_WT_F5_Rep*bam wt-f5-reps 2 SL_WT_F5_Rep*bam wt-f5-reps
GOMAXPROCS=6 parallel -j1 -N4 overlap-end -f {1} -pair {2} -pair {3} -out overlap-data/{4} ::: 0 SL_WT_F2_Rep*bam wt-f2-reps 1 SL_WT_F2_Rep*bam wt-f2-reps 2 SL_WT_F2_Rep*bam wt-f2-reps

# generate data for end overlap analysis denesting long
GOMAXPROCS=6 parallel -j1 -N5 overlap-end -denest -f {1} -pair {2} -pair {3} -pair {4} -out overlap-data-denest/{5} ::: 0 SL_Mut_F2_Rep*bam mut-f2-reps 1 SL_Mut_F2_Rep*bam mut-f2-reps 2 SL_Mut_F2_Rep*bam mut-f2-reps 0 SL_Mut_F5_Rep*bam mut-f5-reps 1 SL_Mut_F5_Rep*bam mut-f5-reps 2 SL_Mut_F5_Rep*bam mut-f5-reps 0 SL_WT_F5_Rep*bam wt-f5-reps 1 SL_WT_F5_Rep*bam wt-f5-reps 2 SL_WT_F5_Rep*bam wt-f5-reps
GOMAXPROCS=6 parallel -j1 -N4 overlap-end -denest -f {1} -pair {2} -pair {3} -out overlap-data-denest/{4} ::: 0 SL_WT_F2_Rep*bam wt-f2-reps 1 SL_WT_F2_Rep*bam wt-f2-reps 2 SL_WT_F2_Rep*bam wt-f2-reps

# generate data for end overlap analysis denesting long and considering only fully containing long
GOMAXPROCS=6 parallel -j1 -N5 overlap-end -denest -contain -f {1} -pair {2} -pair {3} -pair {4} -out overlap-data-denest-contained/{5} ::: 0 SL_Mut_F2_Rep*bam mut-f2-reps 1 SL_Mut_F2_Rep*bam mut-f2-reps 2 SL_Mut_F2_Rep*bam mut-f2-reps 0 SL_Mut_F5_Rep*bam mut-f5-reps 1 SL_Mut_F5_Rep*bam mut-f5-reps 2 SL_Mut_F5_Rep*bam mut-f5-reps 0 SL_WT_F5_Rep*bam wt-f5-reps 1 SL_WT_F5_Rep*bam wt-f5-reps 2 SL_WT_F5_Rep*bam wt-f5-reps
GOMAXPROCS=6 parallel -j1 -N4 overlap-end -denest -contain -f {1} -pair {2} -pair {3} -out overlap-data-denest-contained/{4} ::: 0 SL_WT_F2_Rep*bam wt-f2-reps 1 SL_WT_F2_Rep*bam wt-f2-reps 2 SL_WT_F2_Rep*bam wt-f2-reps
