# Upset plot

library(UpSetR)
library(venneuler)
library(eulerr)

setwd("/Users/jessicanordlund/Dropbox/Research/EpiQC/")

my.assay.cols <- c("MethylSeq" = "#377eb8", 
                   "TruSeq" = "#e41a1c",
                   "EMSeq" = "#4daf4a",
                   "SPLAT" = "#984ea3",
                   "EPIC" = "#ff7f00",
                   "Arrays" = "#ffff33",
                   "TrueMethyl" = "#a65628")

HG <- read.table(gzfile("/Users/jessicanordlund/Downloads/EPIQC/upset_HG002_LAB01_10x.txt.gz"), 
                      header= TRUE)   

#Remove EPIC
indx <- grepl('EPIC', colnames(HG))
test <- HG[, !indx]

# Change name
colnames(test)[grepl('OXBS', colnames(test))] <- "TrueMethyl"
colnames(test) <- gsub('\\_.*$','', colnames(test))

#Remove TrueMethyl (bec using 20x it is missing)
indx <-  grepl('TrueMethyl', colnames(test))
test <- test[,!indx]


# Eulerplot
fit1 <- euler(test)
plot(fit1,
     edges = my.assay.cols[match(names(test), names(my.assay.cols))],
     #fills = paste(my.assay.cols[match(names(test), names(my.assay.cols))], "20", sep=""),
     fills= NA,
     fontsize = 8,
     quantities = list(fontsize = 8))


# Upset plot
upset(test, order.by="freq",
      sets.bar.color = my.assay.cols[match(names(test), names(my.assay.cols))], 
      point.size = 3.5, line.size = 1, 
      sets.x.label = "CpG sites", 
      text.scale = c(2, 1.3, 2, 1, 2, 0.75))




