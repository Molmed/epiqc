library(data.table)
library(GenomicRanges)
library(SummarizedExperiment)

args <- commandArgs(trailingOnly=T)

beds <- args[1:(length(args) - 3)]
nsamp <- length(beds)

load(args[length(args) - 2]) #hg38 cg probe_coords

cov <- args[length(args) - 1]

upset <- matrix(0, nrow=length(probe_coords), ncol=nsamp)

for (i in 1:nsamp) {
	fn <- beds[i]
	print(fn)
	tmp <- fread(fn, skip='chr')
	gr <- GRanges(seqnames=tmp$V1, ranges=IRanges(start=tmp$V2, end=tmp$V3))
	print('made granges')
	hits <- findOverlaps(probe_coords, gr, ignore.strand=T)
	print('found hits')
	tmp <- tmp[subjectHits(hits)]
	cc <- tmp$V5 + tmp$V6
	rm(gr, tmp)
	gc()

	upset[queryHits(hits),i] <- as.integer(cc >= as.numeric(cov))
}

prep <- sapply(strsplit(basename(beds), '[_.]'), '[', 2)
sample <- sapply(strsplit(basename(beds), '[_.]'), '[', 3)
lab <- sapply(strsplit(basename(beds), '[_.]'), '[', 4)
colnames(upset) <- paste(prep, sample, lab, cov, sep='_')

fwrite(upset, file=args[length(args)], sep='\t', quote=F, row.names=F, compress='gzip')

