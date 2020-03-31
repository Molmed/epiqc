library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)

args <- commandArgs(trailingOnly=T)

load(args[1])

tmp <- fread(args[2], skip='chr')

gr <- GRanges(seqnames=tmp$V1, ranges=IRanges(start=tmp$V2 + 1, end=tmp$V3))
hits <- findOverlaps(probe_coords, gr, ignore.strand=T)
tmp <- tmp[subjectHits(hits)]
gc()
tmp$site_idx <- queryHits(hits)
dinuc_sum <- tmp[,.(meth=sum(V5), unmeth=sum(V6)), by=site_idx]
dinuc_sum[,':='(chr=as.character(seqnames(probe_coords))[site_idx], 
				start=start(probe_coords)[site_idx], 
				end=end(probe_coords)[site_idx])]
rm(tmp, gr, hits)
gc()
fwrite(dinuc_sum[,.(chr, start, end, beta=meth/(meth+unmeth), meth, unmeth)], 
	file=args[3], quote=F, row.names=F, col.names=F, sep='\t', compress='gzip')
