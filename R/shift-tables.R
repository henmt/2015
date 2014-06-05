require(reshape)

truncShiftTable <- function(wtFile, mutFile, tissue, type) {
	span <- c(26:36)
	offset <- 36

	mut <- read.csv(mutFile)
	mut.scale <- rbind(
		mut[2,span]/sum(mut[2,span]),
		mut[4,span]/sum(mut[4,span]),
		mut[6,span]/sum(mut[6,span])
	)
	mut.scale <- cbind("Mutant", mut.scale)
	names(mut.scale) <- c("Genotype", span-offset)

	cat("\nSummary mut ends\n\n")
	print(summary(mut.scale))

	wt <- read.csv(wtFile)
	wt.scale <- rbind(
		wt[2,span]/sum(wt[2,span]),
		wt[4,span]/sum(wt[4,span]),
		wt[6,span]/sum(wt[6,span])
	)
	wt.scale <- cbind("Wild-Type", wt.scale)
	names(wt.scale) <- c("Genotype", span-offset)

	cat("\nSummary wt ends\n\n")
	print(summary(wt.scale))
}

# 1° piRNA strict:
truncShiftTable("strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv")
truncShiftTable("strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv")

# 2° piRNA strict:
truncShiftTable("strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv")
truncShiftTable("strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv")

# 1° piRNA not-strict:
truncShiftTable("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "not-strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv")
truncShiftTable("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "not-strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv")

# 2° piRNA not-strict:
truncShiftTable("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "not-strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv")
truncShiftTable("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "not-strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv")


# for contained ends only - only valid if overlap-end is run with -contain
endTruncTableContained <- function(file, tissue, type) {
	ends <- read.csv(file)

	span <- c(26:36)
	offset <- 36

	five.scale <- rbind(
		ends[1,span]/sum(ends[1,span]),
		ends[3,span]/sum(ends[3,span]),
		ends[5,span]/sum(ends[5,span])
	)
	five.scale <- cbind("5'", five.scale)
	names(five.scale) <- c("End", span-offset)

	cat("\nSummary 5' end\n\n")
	print(summary(five.scale))

	three.scale <- rbind(
		ends[2,span]/sum(ends[2,span]),
		ends[4,span]/sum(ends[4,span]),
		ends[6,span]/sum(ends[6,span])
	)
	three.scale <- cbind("3'", three.scale)
	names(three.scale) <- c("End", span-offset)

	cat("\nSummary 3' end\n\n")
	print(summary(three.scale))
}

# 1° piRNA strict:
endTruncTableContained("strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv")
endTruncTableContained("strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv")

# 2° piRNA strict:
endTruncTableContained("strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv")
endTruncTableContained("strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv")

# 1° piRNA not-strict:
endTruncTableContained("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv")
endTruncTableContained("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv")

# 2° piRNA not-strict:
endTruncTableContained("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv")
endTruncTableContained("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv")

