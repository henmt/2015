require(reshape)
require(ggplot2)
require(RColorBrewer)

truncShiftPlot <- function(wtFile, mutFile, tissue, type, out) {
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

	wt <- read.csv(wtFile)
	wt.scale <- rbind(
		wt[2,span]/sum(wt[2,span]),
		wt[4,span]/sum(wt[4,span]),
		wt[6,span]/sum(wt[6,span])
	)
	wt.scale <- cbind("Wild-Type", wt.scale)
	names(wt.scale) <- c("Genotype", span-offset)

	all.scale <- melt(rbind(wt.scale, mut.scale))
	all.scale$Genotype <- factor(all.scale$Genotype, levels=c("Wild-Type", "Mutant"))

	p <- ggplot(all.scale, aes(factor(variable), value))
	p <- p + geom_boxplot(aes(fill = Genotype)) + scale_fill_manual(values = rev(brewer.pal(3, "Set1")[c(1:2)])) + ggtitle(paste("3' Truncation", tissue, type)) + theme(plot.title = element_text(size=16, face="bold", vjust=1.5), axis.title = element_text(size=14)) + xlab("Length change at 3' end") + ylab("Proportion of short piRNA features")
	ggsave(out, plot = p)
}

# 1° piRNA:
truncShiftPlot("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "not-strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "Spermatocyte", "1° piRNA", "f2-primary.pdf")
truncShiftPlot("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "not-strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "Spermatid", "1° piRNA", "f5-primary.pdf")

# 2° piRNA:
truncShiftPlot("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "not-strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "Spermatocyte", "2° piRNA", "f2-secondary.pdf")
truncShiftPlot("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "not-strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "Spermatid", "2° piRNA", "f5-secondary.pdf")


endTruncPlot <- function(file, tissue, type, out) {
	ends <- read.csv(file)

	span <- c(5:62)
	offset <- 36

	five.scale <- rbind(
		ends[1,span]/sum(ends[1,span]),
		ends[3,span]/sum(ends[3,span]),
		ends[5,span]/sum(ends[5,span])
	)
	five.scale <- cbind("5'", five.scale)
	names(five.scale) <- c("End", span-offset)

	three.scale <- rbind(
		ends[2,span]/sum(ends[2,span]),
		ends[4,span]/sum(ends[4,span]),
		ends[6,span]/sum(ends[6,span])
	)
	three.scale <- cbind("3'", three.scale)
	names(three.scale) <- c("End", span-offset)

	all.scale <- melt(rbind(five.scale, three.scale))
	all.scale$End <- factor(all.scale$End, levels=c("5'", "3'"))

	p <- ggplot(all.scale, aes(factor(variable), value), size=0.5)
	p <- p + geom_boxplot(aes(fill = End), linetype=0) + scale_fill_manual(values = rev(brewer.pal(3, "Set1")[c(1:2)])) + ggtitle(paste("End Truncation", tissue, type)) + theme(plot.title = element_text(size=16, face="bold", vjust=1.5), axis.title = element_text(size=14)) + xlab("Length change") + ylab("Proportion of short piRNA features") + scale_x_discrete(breaks=c(-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 25))
	ggsave(out, plot = p)
}

# 1° piRNA:
endTruncPlot("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "WT Spermatocyte", "1° piRNA", "ends-f2-primary.pdf")
endTruncPlot("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "WT Spermatid", "1° piRNA", "ends-f5-primary.pdf")
endTruncPlot("not-strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "Mutant Spermatocyte", "1° piRNA", "mut-ends-f2-primary.pdf")
endTruncPlot("not-strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "Mutant Spermatid", "1° piRNA", "mut-ends-f5-primary.pdf")

# 2° piRNA:
endTruncPlot("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "WT Spermatocyte", "2° piRNA", "ends-f2-secondary.pdf")
endTruncPlot("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "WT Spermatid", "2° piRNA", "ends-f5-secondary.pdf")
endTruncPlot("not-strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "Mutant Spermatocyte", "2° piRNA", "mut-ends-f2-secondary.pdf")
endTruncPlot("not-strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "Mutant Spermatid", "2° piRNA", "mut-ends-f5-secondary.pdf")


# for contained ends only - only valid if overlap-end is run with -contain
endTruncPlotContained <- function(file, tissue, type, out) {
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

	three.scale <- rbind(
		ends[2,span]/sum(ends[2,span]),
		ends[4,span]/sum(ends[4,span]),
		ends[6,span]/sum(ends[6,span])
	)
	three.scale <- cbind("3'", three.scale)
	names(three.scale) <- c("End", span-offset)

	all.scale <- melt(rbind(five.scale, three.scale))
	all.scale$End <- factor(all.scale$End, levels=c("5'", "3'"))

	p <- ggplot(all.scale, aes(factor(variable), value), size=0.5)
	p <- p + geom_boxplot(aes(fill = End)) + scale_fill_manual(values = rev(brewer.pal(3, "Set1")[c(1:2)])) + ggtitle(paste("End Truncation", tissue, type)) + theme(plot.title = element_text(size=16, face="bold", vjust=1.5), axis.title = element_text(size=14)) + xlab("Length change") + ylab("Proportion of short piRNA features")
	ggsave(out, plot = p)
}

# 1° piRNA:
endTruncPlotContained("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "WT Spermatocyte", "1° piRNA", "ends-f2-primary-contain.pdf")
endTruncPlotContained("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "WT Spermatid", "1° piRNA", "ends-f5-primary-contain.pdf")
endTruncPlotContained("not-strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "Mutant Spermatocyte", "1° piRNA", "mut-ends-f2-primary-contain.pdf")
endTruncPlotContained("not-strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-U1.csv", "Mutant Spermatid", "1° piRNA", "mut-ends-f5-primary-contain.pdf")

# 2° piRNA:
endTruncPlotContained("not-strict/SL_WT_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "WT Spermatocyte", "2° piRNA", "ends-f2-secondary-contain.pdf")
endTruncPlotContained("not-strict/SL_WT_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_WT_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "WT Spermatid", "2° piRNA", "ends-f5-secondary-contain.pdf")
endTruncPlotContained("not-strict/SL_Mut_F2_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F2_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "Mutant Spermatocyte", "2° piRNA", "mut-ends-f2-secondary-contain.pdf")
endTruncPlotContained("not-strict/SL_Mut_F5_Rep1_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep2_R1.trimmed_HiSeq2000.bwa.mm10-SL_Mut_F5_Rep3_R1.trimmed_HiSeq2000.bwa.mm10-A10.csv", "Mutant Spermatid", "2° piRNA", "mut-ends-f5-secondary-contain.pdf")

