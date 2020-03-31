library(data.table)
library(GenomicRanges)
library(SummarizedExperiment)

args <- commandArgs(trailingOnly=T)

beds <- args[1:(length(args) - 3)]
nsamp <- length(beds)

load(args[length(args) - 2]) #array_anot_with_hg38.Rd

anot <- anot[anot$isCG38,]
probe_coords <- GRanges(seqnames=anot$chr38, ranges=IRanges(start=anot$pos38, width=2), name=anot$Name, chr_hg19=anot$chr, pos_hg19=anot$pos)
rm(anot)
gc()

meth <- unmeth <- matrix(NA, nrow=length(probe_coords), ncol=nsamp)

for (i in 1:nsamp) {
	fn <- beds[i]
	print(fn)
	tmp <- fread(fn, skip='chr')
	gr <- GRanges(seqnames=tmp$V1, ranges=IRanges(start=tmp$V2 + 1, end=tmp$V3))
	print('made granges')
	hits <- findOverlaps(probe_coords, gr, ignore.strand=T)
	print('found hits')
	tmp <- tmp[subjectHits(hits)]
	rm(gr)
	gc()

	meth[queryHits(hits),i] <- tmp$V5
	unmeth[queryHits(hits),i] <- tmp$V6

	rm(tmp)
	gc()
}

prep <- sapply(strsplit(basename(beds), '[_.]'), '[', 2)
sample <- sapply(strsplit(basename(beds), '[_.]'), '[', 3)
lab <- sapply(strsplit(basename(beds), '[_.]'), '[', 4)

cold <- DataFrame(Sample=sample, Prep=prep, Lab=lab, row.names=basename(beds))
bs_se <- SummarizedExperiment(assays=list(Meth=meth, Unmeth=unmeth), rowRanges=probe_coords, colData=cold)

save(bs_se, file=args[length(args) - 1])

mapping <- data.table(chr_hg38=as.character(seqnames(probe_coords)), pos_hg38=start(probe_coords), 
			name=probe_coords$name, chr_hg19=probe_coords$chr_hg19, pos_hg19=probe_coords$pos_hg19)

fwrite(mapping, file=args[length(args)], sep='\t', quote=F, row.names=F)

