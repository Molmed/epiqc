# Upset plot

library(UpSetR)
library(venneuler)
library(eulerr)


my.assay.cols <- c("MethylSeq" = "#377eb8", #1
                   "TruSeq" = "#e41a1c", #2
                   "EMSeq" = "#4daf4a", #3
                   "SPLAT" = "#984ea3",#4
                   "EPIC" = "#ff7f00", #5
                   "Arrays" = "#ffff33",#6
                   "TrueMethyl" = "#a65628")#7

HG <- read.table(gzfile("/Users/jessicanordlund/Downloads/EPIQC/10xDS/upset_HG001_LAB01_5x.txt.gz"), 
                 header= TRUE)   
HG <- read.table(gzfile("/Users/jessicanordlund/Downloads/EPIQC/10xDS/upset_HG002_LAB01_5x.txt.gz"), 
                      header= TRUE)   
HG <- read.table(gzfile("/Users/jessicanordlund/Downloads/EPIQC/10xDS/upset_HG003_LAB01_5x.txt.gz"), 
                 header= TRUE)  
HG <- read.table(gzfile("/Users/jessicanordlund/Downloads/EPIQC/10xDS/upset_HG004_LAB01_5x.txt.gz"), 
                 header= TRUE)  
HG <- read.table(gzfile("/Users/jessicanordlund/Downloads/EPIQC/10xDS/upset_HG005_LAB01_5x.txt.gz"), 
                 header= TRUE)  
HG <- read.table(gzfile("/Users/jessicanordlund/Downloads/EPIQC/10xDS/upset_HG006_LAB01_5x.txt.gz"), 
                 header= TRUE)   
HG <- read.table(gzfile("/Users/jessicanordlund/Downloads/EPIQC/10xDS/upset_HG007_LAB01_5x.txt.gz"), 
                 header= TRUE)  

#Remove EPIC
indx <- grepl('EPIC', colnames(HG))
test <- HG[, !indx]

# Change name
colnames(test)[grepl('OXBS', colnames(test))] <- "TrueMethyl"
colnames(test) <- gsub('\\_.*$','', colnames(test))

HG001.col <- rev(my.assay.cols[c(2,4,1,3,7)])
HG002.col <- rev(my.assay.cols[c(2,4,1,7,3)])
HG003.col <- rev(my.assay.cols[c(2,4,1,7,3)])
HG004.col <- rev(my.assay.cols[c(2,4,1,7,3)])
HG005.col <- rev(my.assay.cols[c(2,4,1,7,3)])
HG006.col <- rev(my.assay.cols[c(2,4,1,7,3)])
HG007.col <- rev(my.assay.cols[c(2,4,1,7,3)])


# For supplementary
#c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
pdf(file="Upset_HG007_10xDS.pdf", width=9, height=5)
upset(test, 
      sets.bar.color = HG007.col,
      order.by="freq",
      point.size = 2, 
      line.size = 1, 
      number.angles = 0,
      sets.x.label = "N CpG sites", 
      text.scale = c(1.2, 1.2, 1.2, 1.2, 1.1, 0.7))
dev.off()




# Upset plot- for main fig

pdf(file="Upset_HG002_10xDS_main.pdf", width=7, height=6)
upset(test, 
      nintersects = 10,
      sets.bar.color = HG002.col,
      order.by="freq",
      point.size = 3.5, 
      line.size = 1, 
      number.angles = 0,
      sets.x.label = "N CpG sites", 
      text.scale = c(2, 1.4, 2, 1.4, 2, 1.2))

dev.off()



# Eulerplot
fit1 <- euler(test)
plot(fit1,
     edges = my.assay.cols[match(names(test), names(my.assay.cols))],
     #fills = paste(my.assay.cols[match(names(test), names(my.assay.cols))], "20", sep=""),
     fills= NA,
     fontsize = 8,
     quantities = list(fontsize = 8))



