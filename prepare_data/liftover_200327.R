library(minfi)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)

args <- commandArgs(trailingOnly=T)

load(args[1]) # rgset.RData
anot <- getAnnotation(rgset)

anot_gr <- GRanges(seqnames=anot$chr, ranges=IRanges(start=anot$pos, width=1), strand=anot$strand)
names(anot_gr) <- rownames(anot)

ch <- import.chain(args[2]) # hg19ToHg38.over.chain

anot_gr38 <- unlist(liftOver(anot_gr, ch))

strand38 <- pos38 <- chr38 <- rep(NA, length(anot))

idx <-  match(rownames(anot), names(anot_gr38))

chr38[!is.na(idx)] <- seqnames(anot_gr38)[idx[!is.na(idx)]]
pos38[!is.na(idx)] <- start(anot_gr38)[idx[!is.na(idx)]]
strand38[!is.na(idx)] <- strand(anot_gr38)[idx[!is.na(idx)]]

ds <- which(anot$strand != strand38)
pos38[ds] <- pos38[ds] - 1

anot$chr38 <- chr38
anot$pos38 <- pos38
anot$strand38 <- strand38

was_lifted <- !is.na(anot$pos38)

tmp_gr <- resize(GRanges(seqnames=anot$chr38[was_lifted], ranges=IRanges(start=anot$pos38[was_lifted], width=2), strand='+'), 6, fix='center')
ss <- as.character(getSeq(Hsapiens, tmp_gr))
ssCG <- lapply(gregexpr('CG', ss), function(x) as.vector(x) - 3)

anot$off38 <- NA
anot$off38[was_lifted] <- unstrsplit(CharacterList(ssCG), ',')
anot$isCG38 <- F
anot$isCG38[grep('0', anot$off38)] <- T

tt <- table(offset=anot$off38, strand=anot$strand, strand38=anot$strand38)

save(anot, file=args[3]) # array_anot_with_hg38.Rd
