# EPIQC CPG sites by Coverage plots
# Jessica Nordlund - March 2020

#Load libs
library(tidyr)
library(ggplot2)

#Set directors and grab downsampled coverage files (from file_coverage.R (see prepare_data)
setwd("/Users/jessicanordlund/Dropbox/Research/EpiQC/")
test.data <- read.csv("cov_200215_DS.txt", header= T, sep="\t") 

# Format data structure
dat <- gather(test.data[,c(1:9)], variable, value, -Prep)
dat <- dat[-grep("HG", dat$value),]
dat$value <- as.numeric(dat$value)
dat$value <- round(dat$value/1000000, digits= 2)
dat$Library_type <- factor(dat$Prep,levels = c("EMSeq", "MethylSeq", "SPLAT", "TruSeq", "EPIC"))
dat$variable <-  factor(dat$variable, levels = c("cov_1x", "cov_5x", "cov_10x", "cov_20x", "cov_30x", "cov_50x", "cov_100x")) 


#DOT-BOX PLOT 
ggplot(dat, aes(x = variable, y = value, fill = Library_type)) +
  geom_boxplot(outlier.size = 0, alpha= 0.5) +
  geom_point(pch = 21, position = position_jitterdodge(), size = 3, alpha=0.5) +
  xlab("minimum read coverage") +
  ylab("million CpG sites detected")
