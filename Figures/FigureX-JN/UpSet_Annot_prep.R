# module load R/3.6.1
# module load R_packages/3.6.1

library(data.table)
library(dplyr)
library(UpSetR)

setwd('/proj/uppstore2017167/EPIQC/Homer_annotations')
anot <- fread('/proj/uppstore2017167/EPIQC/upset_coords.bed.gz')
anot2 <- fread('/proj/uppstore2017167/EPIQC/cg_annotation/anot/CG_sites_100bp_anot_table.gz')
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

#------------------------------------------------#
# REMOVE X AND Y CHR
#------------------------------------------------#

my.emseq <- my.emseq[with(anot, !V1 %in% c("chrY", "chrX")) ,]
my.methylseq <- my.methylseq[with(anot, !V1 %in% c("chrY", "chrX")) ,]
my.truseq <- my.truseq[with(anot, !V1 %in% c("chrY", "chrX")) ,]
my.truemethyl <- my.truemethyl[with(anot, !V1 %in% c("chrY", "chrX")) ,]
my.splat <- my.splat[with(anot, !V1 %in% c("chrY", "chrX")) ,]

anot2 <- anot2[with(anot, !V1 %in% c("chrY", "chrX")) ,]
anot <- anot[with(anot, !V1 %in% c("chrY", "chrX")) ,]
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

# Write files for CPGs "missed" at 5X cov in all libraries aka "MISSED"
# and for CpGs "catched" at 5X or greater in most libraries aka "Catched"
my.emseq[, sums := rowSums(.SD)]
my.splat[, sums := rowSums(.SD)]
my.methylseq[, sums := rowSums(.SD)]
my.truemethyl[, sums := rowSums(.SD)]
my.truseq[, sums := rowSums(.SD)]

fwrite(anot[with(my.emseq, sums < 1),], file="emseq_coords_missed.bed", sep="\t", col.names=F)
fwrite(anot[with(my.emseq, sums >= 4),], file="emseq_coords_catched.bed", sep="\t", col.names=F)

fwrite(anot[with(my.splat, sums < 1),], file="splat_coords_missed.bed", sep="\t",col.names=F)
fwrite(anot[with(my.splat, sums >= 4),], file="splat_coords_catched.bed", sep="\t",col.names=F)

fwrite(anot[with(my.methylseq, sums < 1),], file="methylseq_coords_missed.bed", sep="\t", col.names=F)
fwrite(anot[with(my.methylseq, sums >= 4),], file="methylseq_coords_catched.bed", sep="\t", col.names=F)

fwrite(anot[with(my.truemethyl, sums < 1),], file="truemethyl_coords_missed.bed", sep="\t",col.names=F)
fwrite(anot[with(my.truemethyl, sums >= 4),], file="truemethyl_coords_catched.bed", sep="\t",col.names=F)

fwrite(anot[with(my.truseq, sums < 1),], file="truseq_coords_missed.bed", sep="\t",col.names=F)
fwrite(anot[with(my.truseq, sums >= 4),], file="truseq_coords_catched.bed", sep="\t",col.names=F)

# EM-Seq 
#missed:3254409
#catched:15908511

# MethylSeq
#missed:3969935
#catched:15741190

# Splat
#missed:5542927
#catched:14139319

#Truemethyl
#missed:3124316
#catched: 15801320

#TruSeq
#missed:10157283
#catched:12063788



tmp <- anot2[with(my.truseq, sums < 1),]






# NEW Upsets for: Never caught, usually caught and always caught
merged.dat <- data.table(emseq = select(my.emseq, sums), 
           methylseq= select(my.methylseq, sums),
           splat= select(my.splat, sums),
           truemethyl=select(my.truemethyl, sums), 
           truseq=select(my.truseq, sums))

# Never caught by any rep (genome)
dt.never <- merged.dat
dt.never[dt.never == 0] <- 8
dt.never[dt.never <= 7] <- 0
dt.never[dt.never == 8] <- 1

# Always caught by all reps (genomes)
dt.always <- merged.dat
dt.always[dt.always < 7] <- 0
dt.always[dt.always == 7] <- 1

# Usually caught by most reps 3/7 or more (genomes)
dt.usually <- merged.dat
dt.usually[dt.usually <= 2] <- 0
dt.usually[dt.usually > 2] <- 1


# Save some upset plots of these 
pdf(file="Upset_NeverCaught_10xDS.pdf", width=9, height=6)
upset(dt.never, 
      nsets =  5,
      order.by="freq",
      point.size = 2, 
      line.size = 1, 
      number.angles = 0,
      sets.x.label = "N CpG sites", 
      text.scale = c(1.2, 1.2, 1.2, 1.2, 1.1, 0.7))
dev.off()

pdf(file="Upset_AlwaysCaught_10xDS.pdf", width=9, height=6)
upset(dt.always, 
      nsets =  5,
      order.by="freq",
      point.size = 2, 
      line.size = 1, 
      number.angles = 0,
      sets.x.label = "N CpG sites", 
      text.scale = c(1.2, 1.2, 1.2, 1.2, 1.1, 0.7))
dev.off()

pdf(file="Upset_UsuallyCaught_10xDS.pdf", width=9, height=6)
upset(dt.usually, 
      nsets =  5,
      order.by="freq",
      point.size = 2, 
      line.size = 1, 
      number.angles = 0,
      sets.x.label = "N CpG sites", 
      text.scale = c(1.2, 1.2, 1.2, 1.2, 1.1, 0.7))
dev.off()


# Save annotations
#1 always missed by TruSeq compared to caught by the others ! OBS THIS NOT DONE
fwrite(anot[with(merged.dat, sums < 1),], file="TruSeq_specific_coords.bed", sep="\t", col.names=F)
fwrite(anot[with(my.emseq, sums >= 4),], file="emseq_coords_catched.bed", sep="\t", col.names=F)



# 2 always missed by Splat compared to the others:


