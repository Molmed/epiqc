# Merge replicates, using total count of cytosines across dinucleotide as coverage measure

library(data.table)


args <- commandArgs(trailingOnly=T)

flist <- args[1:(length(args)-1)]
out <- args[length(args)]

print(paste('Merging replicate strandmerged files:\n', paste(flist, collapse='\n')))

print(paste('reading', flist[1]))
dt <- fread(flist[1], drop='V4', skip='chr')
print(paste('Number of sites:', nrow(dt)))
print(paste('Meth Cs:', sum(dt$V5), ', unmeth Cs:', sum(dt$V6)))

if (length(flist) > 1) {
	for (j in 2:length(flist)) {
		print(paste('reading', flist[j]))
		dt_tmp <- fread(flist[j], drop='V4', skip='chr')
		print(paste('Number of sites:', nrow(dt_tmp)))
		print(paste('Meth Cs:', sum(dt_tmp$V5), ', unmeth Cs:', sum(dt_tmp$V6)))
		dt <- merge(dt, dt_tmp, by=c('V1', 'V2', 'V3'), all=T)
		rm(dt_tmp)
		gc()
		dt[,V5:=rowSums(.SD, na.rm=T), .SDcols=grep('^V5\\.', names(dt))]
		dt[,V6:=rowSums(.SD, na.rm=T), .SDcols=grep('^V6\\.', names(dt))]
		dt <- dt[,paste0('V', c(1:3, 5,6))]
		print(paste('Number of merged sites:', nrow(dt)))
		print(paste('Meth Cs:', sum(dt$V5), ', unmeth Cs:', sum(dt$V6)))
		gc()
	}
}
dt[,V4:=V5/(V5+V6)*100]
dt <- dt[,paste0('V', 1:6)]

fwrite(dt, file=out, compress='gzip', sep='\t', quote=F, row.names=F, col.names=F)
print(paste('Wrote replicate merged file:', out))
