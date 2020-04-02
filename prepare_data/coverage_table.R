library(data.table)


args <- commandArgs(trailingOnly=T)

file_in <- args[1]
file_out <- args[2]

bn <- basename(file_in)
print(bn)
tech <- sapply(strsplit(bn, '_'), '[', 2)
sml <- sapply(strsplit(bn, '[_]'), function(x) x[3] )
lab <- sapply(strsplit(bn, '[_]'), function(x) x[4] )
ds <- sapply(strsplit(bn, '[_.]'), function(x) x[5] )

dd <- fread(file_in, skip='chr')

if (nrow(dd) > 0) {
	dd[,cov:=(V5+V6)]
	brk <- c(1, 2, 5, 10, 20, 30, 40, 50, 100, 200, Inf)
	dd[, bins:=cut(cov, brk, right=F, ordered_result=T)]

	foo <- dd[,.(count=.N), by=bins]
	setorder(foo, bins)
	foo[,cumulative:=rev(cumsum(rev(count)))]
	dt <- cbind(data.table(prep=tech, sample=sml, lab=lab, downsampled=ds, source=bn), foo)
	write.table(dt, file=file_out, sep='\t', quote=F, row.names=F, col.names=F)
} else {
	system(paste('touch', file_out))
}
