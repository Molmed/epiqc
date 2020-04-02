library(data.table)

args <- commandArgs(trailingOnly=T)

file_in <- args[1]
file_out <- args[2]

bn <- basename(file_in)
print(bn)
tech <- sapply(strsplit(bn, '_'), '[', 2)
sml <- sapply(strsplit(bn, '_'), function(x) x[3] )
lab <- sapply(strsplit(bn, '[_.]'), function(x) x[4] )

dd <- fread(file_in, skip='chr')


if (nrow(dd) > 0) {
	tot_cov <- data.table(prep=tech, sample=sml, lab=lab, totc=dd[,sum(V5+V6)], avg_cov=dd[,mean(V5+V6, na.rm=T)], med_cov=dd[,median((V5+V6), na.rm=T)], source=bn)
	write.table(tot_cov, file=file_out, sep='\t', quote=F, row.names=F, col.names=F)
} else {
	system(paste('touch', file_out))
}
