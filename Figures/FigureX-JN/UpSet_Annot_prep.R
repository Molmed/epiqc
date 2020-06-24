# module load R/3.6.1
# module load R_packages/3.6.1

library(data.table)
library(dplyr)
library(UpSetR)

anot <- fread('/proj/uppstore2017167/EPIQC/upset_coords.bed.gz')

dat.dir <- '/proj/uppstore2017167/EPIQC/Homer_annotations'
my.files <- list.files(dat.dir, recursive=F, full.names= T, pattern= "LAB01_10x.txt.gz")

my.emseq <- data.table()
my.methylseq <-   data.table()
my.truseq <- data.table()
my.splat <- data.table()
my.truemethyl <- data.table()

for (f in my.files) {
  df <- fread(f)
  my.emseq <- cbind(my.emseq, df %>% select(contains("EMSeq")))
  my.methylseq <- cbind(my.methylseq, df %>% select(contains("MethylSeq")))
  my.truseq <- cbind(my.truseq, df %>% select(contains("TruSeq")))
  my.truemethyl <- cbind(my.truemethyl, df %>% select(contains("OXBS")))
  my.splat <- cbind(my.splat, df %>% select(contains("SPLAT")))
}

# EMSeq
my.emseq[, sums := rowSums(.SD)]
fwrite(anot[with(my.emseq, sums >=3),], file="emseq_coords.bed", sep="\t", col.names=F)

my.splat[, sums := rowSums(.SD)]
fwrite(anot[with(my.splat, sums >=3),], file="splat_coords.bed", sep="\t",col.names=F)


#------------------------------------------------#
# UPSET PLOTS
#------------------------------------------------#


pdf(file="Upset_EMSeq_10xDS.pdf", width=9, height=6)
upset(my.emseq, 
      order.by="freq",
      nsets =  7,
      point.size = 2, 
      line.size = 1, 
      number.angles = 0,
      sets.x.label = "N CpG sites", 
      text.scale = c(1.2, 1.2, 1.2, 1.2, 1.1, 0.7))
dev.off()

pdf(file="Upset_SPLAT_10xDS.pdf", width=9, height=6)
upset(my.splat, 
      nsets =  7,
      order.by="freq",
      point.size = 2, 
      line.size = 1, 
      number.angles = 0,
      sets.x.label = "N CpG sites", 
      text.scale = c(1.2, 1.2, 1.2, 1.2, 1.1, 0.7))
dev.off()

pdf(file="Upset_MethylSeq_10xDS.pdf", width=9, height=6)
upset(my.methylseq, 
      nsets =  7,
      order.by="freq",
      point.size = 2, 
      line.size = 1, 
      number.angles = 0,
      sets.x.label = "N CpG sites", 
      text.scale = c(1.2, 1.2, 1.2, 1.2, 1.1, 0.7))
dev.off()

pdf(file="Upset_TruSeq_10xDS.pdf", width=9, height=6)
upset(my.truseq, 
      nsets =  7,
      order.by="freq",
      point.size = 2, 
      line.size = 1, 
      number.angles = 0,
      sets.x.label = "N CpG sites", 
      text.scale = c(1.2, 1.2, 1.2, 1.2, 1.1, 0.7))
dev.off()

pdf(file="Upset_TrueMethyl_10xDS.pdf", width=9, height=6)
upset(my.truemethyl, 
      nsets =  7,
      order.by="freq",
      point.size = 2, 
      line.size = 1, 
      number.angles = 0,
      sets.x.label = "N CpG sites", 
      text.scale = c(1.2, 1.2, 1.2, 1.2, 1.1, 0.7))
dev.off()
