# Strand and replicate merging, bsseq liftover to epic array and downsampling of bsseq bedgraphs

library(BSgenome.Hsapiens.UCSC.hg38)

args <- commandArgs(trailingOnly=T)
output <- args[1]

chrs <- paste0('chr', c(1:22, 'X', 'Y', 'M'))
probe_coords <- unlist(GRangesList(lapply(chrs, function(x) GRanges(seqnames=x, ranges=IRanges(start=start(matchPattern("CG", Hsapiens[[x]])), width=2)))))

save(probe_coords, file=output)
