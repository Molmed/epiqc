library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)

args <- commandArgs(trailingOnly=T)

load(args[1])

tmp <- fread(args[2], skip='chr')

print(paste('Strandmerging', args[2]))
print(paste('Number of cg sites in reference chrs 1-22, X, Y, M:', length(probe_coords)))

gr <- GRanges(seqnames=tmp$V1, ranges=IRanges(start=tmp$V2 + 1, end=tmp$V3))
print(paste('Number of Cs assayed in bedgraph:', length(gr)))
print(paste('Number of Cs assayed in bedgraph chr1-22,X,Y,M:', 
	length(keepSeqlevels(gr, seqlevels(probe_coords),  pruning.mode='coarse'))))
 
hits <- findOverlaps(probe_coords, gr, ignore.strand=T)
print(paste('Number of reference cg sites covered by bedgraph data:', length(unique(queryHits(hits)))))
print(paste('Number of bedgraph Cs intersecting reference cg sites:', length(unique(subjectHits(hits)))))

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

print(paste('Wrote strandmerged file with', nrow(dinuc_sum), 'cg sites'))
print(paste('Methylated Cs:', dinuc_sum[,sum(meth)], ', unmethylated Cs:', dinuc_sum[,sum(unmeth)]))
