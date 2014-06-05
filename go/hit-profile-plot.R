# GOMAXPROCS=4 parallel -N2 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2} > wt-profile-f2.scov' ::: SL_WT_F2*bam
# GOMAXPROCS=4 parallel -N2 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2} -f 1 > wt-profile-f2-U1.scov' ::: SL_WT_F2*bam
# GOMAXPROCS=4 parallel -N2 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2} -f 2 > wt-profile-f2-A10.scov' ::: SL_WT_F2*bam
# GOMAXPROCS=4 parallel -N2 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2} -f 1 -strict > wt-profile-f2-U1-only.scov' ::: SL_WT_F2*bam
# GOMAXPROCS=4 parallel -N2 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2} -f 2 -strict > wt-profile-f2-A10-only.scov' ::: SL_WT_F2*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} > mut-profile-f2.scov' ::: SL_Mut_F2*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 1 > mut-profile-f2-U1.scov' ::: SL_Mut_F2*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 2 > mut-profile-f2-A10.scov' ::: SL_Mut_F2*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 1 -strict > mut-profile-f2-U1-only.scov' ::: SL_Mut_F2*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 2 -strict > mut-profile-f2-A10-only.scov' ::: SL_Mut_F2*bam

# ## requires change of grep pattern for each line type
# grep L1Md_T wt-profile-f2.scov > wt-profile-f2-T.scov
# grep L1Md_T mut-profile-f2.scov > mut-profile-f2-T.scov
# grep L1Md_T wt-profile-f2-U1.scov > wt-profile-f2-U1-T.scov
# grep L1Md_T mut-profile-f2-U1.scov > mut-profile-f2-U1-T.scov
# grep L1Md_T wt-profile-f2-A10.scov > wt-profile-f2-A10-T.scov
# grep L1Md_T mut-profile-f2-A10.scov > mut-profile-f2-A10-T.scov
# grep L1Md_T wt-profile-f2-U1-only.scov > wt-profile-f2-U1-only-T.scov
# grep L1Md_T mut-profile-f2-U1-only.scov > mut-profile-f2-U1-only-T.scov
# grep L1Md_T wt-profile-f2-A10-only.scov > wt-profile-f2-A10-only-T.scov
# grep L1Md_T mut-profile-f2-A10-only.scov > mut-profile-f2-A10-only-T.scov

# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} > wt-profile-f5.scov' ::: SL_WT_F5*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 1 > wt-profile-f5-U1.scov' ::: SL_WT_F5*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 2 > wt-profile-f5-A10.scov' ::: SL_WT_F5*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 1 -strict > wt-profile-f5-U1-only.scov' ::: SL_WT_F5*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 2 -strict > wt-profile-f5-A10-only.scov' ::: SL_WT_F5*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} > mut-profile-f5.scov' ::: SL_Mut_F5*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 1 > mut-profile-f5-U1.scov' ::: SL_Mut_F5*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 2 > mut-profile-f5-A10.scov' ::: SL_Mut_F5*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 1 -strict > mut-profile-f5-U1-only.scov' ::: SL_Mut_F5*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 2 -strict > mut-profile-f5-A10-only.scov' ::: SL_Mut_F5*bam

# ## requires change of grep pattern for each line type
# grep L1Md_A wt-profile-f5.scov > wt-profile-f5-A.scov
# grep L1Md_A mut-profile-f5.scov > mut-profile-f5-A.scov
# grep L1Md_A wt-profile-f5-U1.scov > wt-profile-f5-U1-A.scov
# grep L1Md_A mut-profile-f5-U1.scov > mut-profile-f5-U1-A.scov
# grep L1Md_A wt-profile-f5-A10.scov > wt-profile-f5-A10-A.scov
# grep L1Md_A mut-profile-f5-A10.scov > mut-profile-f5-A10-A.scov
# grep L1Md_A wt-profile-f5-U1-only.scov > wt-profile-f5-U1-only-A.scov
# grep L1Md_A mut-profile-f5-U1-only.scov > mut-profile-f5-U1-only-A.scov
# grep L1Md_A wt-profile-f5-A10-only.scov > wt-profile-f5-A10-only-A.scov
# grep L1Md_A mut-profile-f5-A10-only.scov > mut-profile-f5-A10-only-A.scov


# line coverage for fragments that have at least one piRNA hit
# GOMAXPROCS=4 parallel -N2 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/repeat-hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2} > wt-line-profile-f2.cov' ::: SL_WT_F2*bam
# GOMAXPROCS=4 parallel -N2 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/repeat-hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2} -f 1 > wt-line-profile-f2-U1.cov' ::: SL_WT_F2*bam
# GOMAXPROCS=4 parallel -N2 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/repeat-hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2} -f 2 > wt-line-profile-f2-A10.cov' ::: SL_WT_F2*bam
# GOMAXPROCS=4 parallel -N2 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/repeat-hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2} -f 1 -strict > wt-line-profile-f2-U1-only.cov' ::: SL_WT_F2*bam
# GOMAXPROCS=4 parallel -N2 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/repeat-hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2} -f 2 -strict > wt-line-profile-f2-A10-only.cov' ::: SL_WT_F2*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/repeat-hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} > mut-line-profile-f2.cov' ::: SL_Mut_F2*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/repeat-hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 1 > mut-line-profile-f2-U1.cov' ::: SL_Mut_F2*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/repeat-hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 2 > mut-line-profile-f2-A10.cov' ::: SL_Mut_F2*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/repeat-hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 1 -strict > mut-line-profile-f2-U1-only.cov' ::: SL_Mut_F2*bam
# GOMAXPROCS=4 parallel -N3 '~/Development/src/bitbucket.org/sacgf/lim_2014/go/repeat-hit-profile -annot mm10.feat.gff -class repeat/LINE -reads {1},{2},{3} -f 2 -strict > mut-line-profile-f2-A10-only.cov' ::: SL_Mut_F2*bam

require(ggplot2)

profile <- function(wf, mf, lf, main, out) {
	wt <- cbind("wild-type", read.delim(wf, header=F))
	mt <- cbind("mutant", read.delim(mf, header=F))
	l <- cbind("cover", read.delim(lf, header=F))
	n <- c("type", "class", "name", "position", "weight")
	names(wt) <- n
	names(mt) <- n
	names(l) <- n
	wt$weight <- wt$weight/sum(wt$weight)
	mt$weight <- mt$weight/sum(mt$weight)
	l$weight <- l$weight/sum(l$weight)
	all <- rbind(wt, mt, l)

	m <- ggplot(all, aes(x=position, weight=weight, color=type, group=type))
	m <- m + geom_density(adjust=0.1) + ggtitle(main) + theme(plot.title = element_text(size=18, face="bold", vjust=3), axis.title = element_text(size=16)) + xlab("Position (bp)")  + scale_x_continuous(breaks=c(0:100)*1000)

	ggsave(out, plot = m)
}

profile("wt-profile-f2-A.cov", "mut-profile-f2-A.cov", "profile-line-L1Md_A.cov", "LINE L1 A short RNA mapping profiles", "profile-f2-A.pdf")
profile("wt-profile-f2-Gf.cov", "mut-profile-f2-Gf.cov", "profile-line-L1Md_Gf.cov", "LINE L1 Gf short RNA mapping profiles", "profile-f2-Gf.pdf")
profile("wt-profile-f2-T.cov", "mut-profile-f2-T.cov", "profile-line-L1Md_T.cov", "LINE L1 T short RNA mapping profiles", "profile-f2-T.pdf")
profile("wt-profile-f2.cov", "mut-profile-f2.cov", "profile-line-L1.cov", "LINE L1 short RNA mapping profiles", "profile-f2.pdf")
profile("wt-profile-f2-E.cov", "mut-profile-f2-E.cov", "profile-line-L1ME.cov", "LINE L1 E short RNA mapping profiles", "profile-f2-E.pdf")
profile("wt-profile-f2-3.cov", "mut-profile-f2-3.cov", "profile-line-L1M3.cov", "LINE L1 3 short RNA mapping profiles", "profile-f2-3.pdf")

profile("wt-profile-f2-U1-A.cov", "mut-profile-f2-U1-A.cov", "profile-line-L1Md_A.cov", "LINE L1 A 1° piRNA mapping profiles", "profile-f2-U1-A.pdf")
profile("wt-profile-f2-U1-Gf.cov", "mut-profile-f2-U1-Gf.cov", "profile-line-L1Md_Gf.cov", "LINE L1 Gf 1° piRNA mapping profiles", "profile-f2-U1-Gf.pdf")
profile("wt-profile-f2-U1-T.cov", "mut-profile-f2-U1-T.cov", "profile-line-L1Md_T.cov", "LINE L1 T 1° piRNA mapping profiles", "profile-f2-U1-T.pdf")
profile("wt-profile-f2-U1.cov", "mut-profile-f2-U1.cov", "profile-line-L1.cov", "LINE L1 1° piRNA mapping profiles", "profile-f2-U1.pdf")
profile("wt-profile-f2-U1-E.cov", "mut-profile-f2-U1-E.cov", "profile-line-L1ME.cov", "LINE L1 E 1° piRNA mapping profiles", "profile-f2-U1-E.pdf")
profile("wt-profile-f2-U1-3.cov", "mut-profile-f2-U1-3.cov", "profile-line-L1M3.cov", "LINE L1 3 1° piRNA mapping profiles", "profile-f2-U1-3.pdf")

profile("wt-profile-f2-A10-A.cov", "mut-profile-f2-A10-A.cov", "profile-line-L1Md_A.cov", "LINE L1 A 2° piRNA mapping profiles", "profile-f2-A10-A.pdf")
profile("wt-profile-f2-A10-Gf.cov", "mut-profile-f2-A10-Gf.cov", "profile-line-L1Md_Gf.cov", "LINE L1 Gf 2° piRNA mapping profiles", "profile-f2-A10-Gf.pdf")
profile("wt-profile-f2-A10-T.cov", "mut-profile-f2-A10-T.cov", "profile-line-L1Md_T.cov", "LINE L1 T 2° piRNA mapping profiles", "profile-f2-A10-T.pdf")
profile("wt-profile-f2-A10.cov", "mut-profile-f2-A10.cov", "profile-line-L1.cov", "LINE L1 2° piRNA mapping profiles", "profile-f2-A10.pdf")
profile("wt-profile-f2-A10-E.cov", "mut-profile-f2-A10-E.cov", "profile-line-L1ME.cov", "LINE L1 E 2° piRNA mapping profiles", "profile-f2-A10-E.pdf")
profile("wt-profile-f2-A10-3.cov", "mut-profile-f2-A10-3.cov", "profile-line-L1M3.cov", "LINE L1 3 2° piRNA mapping profiles", "profile-f2-A10-3.pdf")

profile("wt-profile-f2-U1-only-A.cov", "mut-profile-f2-U1-only-A.cov", "profile-line-L1Md_A.cov", "LINE L1 A strictly 1° piRNA mapping profiles", "profile-f2-U1-only-A.pdf")
profile("wt-profile-f2-U1-only-Gf.cov", "mut-profile-f2-U1-only-Gf.cov", "profile-line-L1Md_Gf.cov", "LINE L1 Gf strictly 1° piRNA mapping profiles", "profile-f2-U1-only-Gf.pdf")
profile("wt-profile-f2-U1-only-T.cov", "mut-profile-f2-U1-only-T.cov", "profile-line-L1Md_T.cov", "LINE L1 T strictly 1° piRNA mapping profiles", "profile-f2-U1-only-T.pdf")
profile("wt-profile-f2-U1-only.cov", "mut-profile-f2-U1-only.cov", "profile-line-L1.cov", "LINE L1 strictly 1° piRNA mapping profiles", "profile-f2-U1-only.pdf")
profile("wt-profile-f2-U1-only-E.cov", "mut-profile-f2-U1-only-E.cov", "profile-line-L1ME.cov", "LINE L1 E strictly 1° piRNA mapping profiles", "profile-f2-U1-only-E.pdf")
profile("wt-profile-f2-U1-only-3.cov", "mut-profile-f2-U1-only-3.cov", "profile-line-L1M3.cov", "LINE L1 3 strictly 1° piRNA mapping profiles", "profile-f2-U1-only-3.pdf")

profile("wt-profile-f2-A10-only-A.cov", "mut-profile-f2-A10-only-A.cov", "profile-line-L1Md_A.cov", "LINE L1 A strictly 2° piRNA mapping profiles", "profile-f2-A10-only-A.pdf")
profile("wt-profile-f2-A10-only-Gf.cov", "mut-profile-f2-A10-only-Gf.cov", "profile-line-L1Md_Gf.cov", "LINE L1 Gf strictly 2° piRNA mapping profiles", "profile-f2-A10-only-Gf.pdf")
profile("wt-profile-f2-A10-only-T.cov", "mut-profile-f2-A10-only-T.cov", "profile-line-L1Md_T.cov", "LINE L1 T strictly 2° piRNA mapping profiles", "profile-f2-A10-only-T.pdf")
profile("wt-profile-f2-A10-only.cov", "mut-profile-f2-A10-only.cov", "profile-line-L1.cov", "LINE L1 strictly 2° piRNA mapping profiles", "profile-f2-A10-only.pdf")
profile("wt-profile-f2-A10-only-E.cov", "mut-profile-f2-A10-only-E.cov", "profile-line-L1ME.cov", "LINE L1 E strictly 2° piRNA mapping profiles", "profile-f2-A10-only-E.pdf")
profile("wt-profile-f2-A10-only-3.cov", "mut-profile-f2-A10-only-3.cov", "profile-line-L1M3.cov", "LINE L1 3 strictly 2° piRNA mapping profiles", "profile-f2-A10-only-3.pdf")


## stranded
sprofile <- function(wf, mf, lf, main, out) {
	wt <- cbind("wild-type", read.delim(wf, header=F))
	mt <- cbind("mutant", read.delim(mf, header=F))
	l <- cbind("cover", read.delim(lf, header=F))
	l <- cbind(l, "plus")
	pn <- c("type", "class", "name", "strand", "position", "weight")
	ln <- c("type", "class", "name", "position", "weight", "strand")
	names(wt) <- pn
	names(mt) <- pn
	names(l) <- ln
	wt$weight <- wt$weight/sum(wt$weight)
	mt$weight <- mt$weight/sum(mt$weight)
	l$weight <- l$weight/sum(l$weight)
	all <- rbind(wt, mt, l)
	all$strandtype <- sprintf("%s/%s", all$type, all$strand)


	m <- ggplot(all, aes(x=position, weight=weight, color=strandtype, group=strandtype))
	m <- m + geom_density(adjust=0.1) + ggtitle(main) + theme(plot.title = element_text(size=18, face="bold", vjust=3), axis.title = element_text(size=16)) + xlab("Position (bp)")  + scale_x_continuous(breaks=c(0:100)*1000)

	ggsave(out, plot = m)
}

# f2
sprofile("wt-profile-f2-A.scov", "mut-profile-f2-A.scov", "profile-line-L1Md_A.cov", "LINE L1 A short RNA mapping profiles", "profile-f2-A-stranded.pdf")
sprofile("wt-profile-f2-T.scov", "mut-profile-f2-T.scov", "profile-line-L1Md_T.cov", "LINE L1 T short RNA mapping profiles", "profile-f2-T-stranded.pdf")

sprofile("wt-profile-f2-U1-A.scov", "mut-profile-f2-U1-A.scov", "profile-line-L1Md_A.cov", "LINE L1 A 1° piRNA mapping profiles", "profile-f2-U1-A-stranded.pdf")
sprofile("wt-profile-f2-U1-T.scov", "mut-profile-f2-U1-T.scov", "profile-line-L1Md_T.cov", "LINE L1 T 1° piRNA mapping profiles", "profile-f2-U1-T-stranded.pdf")

sprofile("wt-profile-f2-A10-A.scov", "mut-profile-f2-A10-A.scov", "profile-line-L1Md_A.cov", "LINE L1 A 2° piRNA mapping profiles", "profile-f2-A10-A-stranded.pdf")
sprofile("wt-profile-f2-A10-T.scov", "mut-profile-f2-A10-T.scov", "profile-line-L1Md_T.cov", "LINE L1 T 2° piRNA mapping profiles", "profile-f2-A10-T-stranded.pdf")

sprofile("wt-profile-f2-U1-only-A.scov", "mut-profile-f2-U1-only-A.scov", "profile-line-L1Md_A.cov", "LINE L1 A strictly 1° piRNA mapping profiles", "profile-f2-U1-only-A-stranded.pdf")
sprofile("wt-profile-f2-U1-only-T.scov", "mut-profile-f2-U1-only-T.scov", "profile-line-L1Md_T.cov", "LINE L1 T strictly 1° piRNA mapping profiles", "profile-f2-U1-only-T-stranded.pdf")

sprofile("wt-profile-f2-A10-only-A.scov", "mut-profile-f2-A10-only-A.scov", "profile-line-L1Md_A.cov", "LINE L1 A strictly 2° piRNA mapping profiles", "profile-f2-A10-only-A-stranded.pdf")
sprofile("wt-profile-f2-A10-only-T.scov", "mut-profile-f2-A10-only-T.scov", "profile-line-L1Md_T.cov", "LINE L1 T strictly 2° piRNA mapping profiles", "profile-f2-A10-only-T-stranded.pdf")

# f5
sprofile("wt-profile-f5-A.scov", "mut-profile-f5-A.scov", "profile-line-L1Md_A.cov", "LINE L1 A short RNA mapping profiles", "profile-f5-A-stranded.pdf")
sprofile("wt-profile-f5-T.scov", "mut-profile-f5-T.scov", "profile-line-L1Md_T.cov", "LINE L1 T short RNA mapping profiles", "profile-f5-T-stranded.pdf")

sprofile("wt-profile-f5-U1-A.scov", "mut-profile-f5-U1-A.scov", "profile-line-L1Md_A.cov", "LINE L1 A 1° piRNA mapping profiles", "profile-f5-U1-A-stranded.pdf")
sprofile("wt-profile-f5-U1-T.scov", "mut-profile-f5-U1-T.scov", "profile-line-L1Md_T.cov", "LINE L1 T 1° piRNA mapping profiles", "profile-f5-U1-T-stranded.pdf")

sprofile("wt-profile-f5-A10-A.scov", "mut-profile-f5-A10-A.scov", "profile-line-L1Md_A.cov", "LINE L1 A 2° piRNA mapping profiles", "profile-f5-A10-A-stranded.pdf")
sprofile("wt-profile-f5-A10-T.scov", "mut-profile-f5-A10-T.scov", "profile-line-L1Md_T.cov", "LINE L1 T 2° piRNA mapping profiles", "profile-f5-A10-T-stranded.pdf")

sprofile("wt-profile-f5-U1-only-A.scov", "mut-profile-f5-U1-only-A.scov", "profile-line-L1Md_A.cov", "LINE L1 A strictly 1° piRNA mapping profiles", "profile-f5-U1-only-A-stranded.pdf")
sprofile("wt-profile-f5-U1-only-T.scov", "mut-profile-f5-U1-only-T.scov", "profile-line-L1Md_T.cov", "LINE L1 T strictly 1° piRNA mapping profiles", "profile-f5-U1-only-T-stranded.pdf")

sprofile("wt-profile-f5-A10-only-A.scov", "mut-profile-f5-A10-only-A.scov", "profile-line-L1Md_A.cov", "LINE L1 A strictly 2° piRNA mapping profiles", "profile-f5-A10-only-A-stranded.pdf")
sprofile("wt-profile-f5-A10-only-T.scov", "mut-profile-f5-A10-only-T.scov", "profile-line-L1Md_T.cov", "LINE L1 T strictly 2° piRNA mapping profiles", "profile-f5-A10-only-T-stranded.pdf")
