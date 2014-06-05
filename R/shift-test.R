require(reshape)
require(vcdExtra)

printf <- function(...) cat(sprintf(...))

# takes t.test or wilcox.test etc. in test param.
testTruncShift <- function(wtFile, mutFile, test) {
	mut <- read.csv(mutFile)
	wt <- read.csv(wtFile)
	mut <- cbind("mut", mut[c(2,4,6),c(26:36)])
	names(mut) <- c("type", -10:0)
	wt <- cbind("wt", wt[c(2,4,6),c(26:36)])
	names(wt) <- c("type", -10:0)
	all <- melt(rbind(wt, mut))
	obs <- expand.dft(all, freq="value")

	printf("\ncompare: %s %s\n\n", wtFile, mutFile)
	test(variable~type, data=obs, conf.int=TRUE)
}


# 1° piRNA strict:
# testTruncShift("strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", t.test)
# testTruncShift("strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", t.test)

# 2° piRNA strict:
# testTruncShift("strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", t.test)
# testTruncShift("strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", t.test)

# 1° piRNA not-strict:
testTruncShift("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "not-strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", t.test)
testTruncShift("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "not-strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", t.test)

# 2° piRNA not-strict:
testTruncShift("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "not-strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", t.test)
testTruncShift("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "not-strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", t.test)
