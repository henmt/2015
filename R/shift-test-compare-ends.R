require(reshape)
require(vcdExtra)

printf <- function(...) cat(sprintf(...))

# takes t.test or wilcox.test etc. in test param.
testTruncShift <- function(file, test) {
	data <- read.csv(file)
	five <- cbind("five", data[c(1,3,5),c(26:36)])
	names(five) <- c("type", -10:0)
	three <- cbind("three", data[c(2,4,6),c(26:36)])
	names(three) <- c("type", -10:0)
	all <- melt(rbind(five, three))
	obs <- expand.dft(all, freq="value")

	printf("\ncompare: %s\n\n", file)
	test(variable~type, data=obs, conf.int=TRUE)
}

# wild-type
# 1° piRNA strict:
# testTruncShift("strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", t.test)
# testTruncShift("strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", t.test)

# 2° piRNA strict:
# testTruncShift("strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", t.test)
# testTruncShift("strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", t.test)

# 1° piRNA not-strict:
testTruncShift("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", t.test)
testTruncShift("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", t.test)

# 2° piRNA not-strict:
testTruncShift("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", t.test)
testTruncShift("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", t.test)

# mutant
# 1° piRNA strict:
# testTruncShift("strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", t.test)
# testTruncShift("strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", t.test)

# 2° piRNA strict:
# testTruncShift("strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", t.test)
# testTruncShift("strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", t.test)

# 1° piRNA not-strict:
testTruncShift("not-strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", t.test)
testTruncShift("not-strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", t.test)

# 2° piRNA not-strict:
testTruncShift("not-strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", t.test)
testTruncShift("not-strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", t.test)
