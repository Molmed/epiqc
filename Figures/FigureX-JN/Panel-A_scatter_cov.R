# EPIQC CPG sites by Coverage plots
# Jessica Nordlund - March 2020

###------------###
library(tidyr)
library(ggplot2)
library(plyr)

###------------###
# Set EPIQC color scheme
###------------###
my.assay.cols <- c("MethylSeq" = "#377eb8", 
               "TruSeq" = "#e41a1c",
               "EMSeq" = "#4daf4a",
               "SPLAT" = "#984ea3",
               "EPIC" = "#ff7f00",
               "Arrays" = "#ffff33",
               "TrueMethyl" = "#a65628")

###------------###
#Set directors and grab downsampled coverage files (from (from file_coverage.R)
###------------###
setwd("/Users/jessicanordlund/Dropbox/Research/EpiQC/")
test.data <- read.csv("cumulative_coverage_table_200403.txt", header= T, sep="\t") #Downsampled coverages


###---------------------###
#     SUBSET DATA 
###---------------------###
# OBS!!! Pick either by Biological Rep or merged high cov data

# Data subset 1, biological replicates
dat <- test.data[!test.data$sample %in% "HG001-7" ,]
dat <- dat[!dat$prep %in% "EPIC" ,]
dat.name <- "byBioRep"

# Data subset 2, biological replicates merged and data downsampled
dat <- test.data[test.data$sample %in% "HG001-7" ,]
dat <- dat[!dat$prep %in% "EPIC" ,]
dat.name <- "BioMergedHC"


###---------------------###
#Make some changes to the names and number formatting for the plots
###---------------------###

# 1 Change the name of OXBS to TrueMethyl & then re-order alphabetically
levels(dat$prep)[levels(dat$prep)=="OXBS"] <- "TrueMethyl"
dat$prep <- factor(as.character(dat$prep))


# 2 Make numbers in Million w 2 digits
dat$value <- round(as.numeric(dat$cumulative_cg_site_count)/1000000, digits= 2)
dat$value2 <- round(as.numeric(dat$cg_site_count)/1000000, digits= 2)


# 3 Change lable names and order levels of factors
new.labels <- gsub("\\,.*", "", gsub("\\[", "\u2265", as.character(levels(dat$coverage_interval))))
dat$variable <-  mapvalues(dat$coverage_interval, from= levels(dat$coverage_interval), to= new.labels)
dat$variable <- factor(dat$variable, levels = unique(dat$variable))
dat$downsampled_cov <- factor(dat$downsampled, levels = unique(dat$downsampled))
dat$variable2 <-  gsub(")", "]", dat$coverage_interval)
dat$variable2 <- factor(dat$variable2, levels = unique(dat$variable2))


###---------------------###
#    PLOTS 
###---------------------###

# Set plot colors by assay
tmp.cols <- my.assay.cols[match(levels(dat$prep), names(my.assay.cols))]

# N of columns depending on which dataset used
if (dat.name %in% "byBioRep") {
    n.cols.plot <- 6
    h.plot <- 2.25
  } else {
    n.cols.plot <- 3
    h.plot <- 4
  }
  
p2 <- ggplot(dat, aes(x = variable, y = value, group = prep)) +
  geom_point((aes(shape=prep, color=prep, size=prep)), 
              position = position_jitter(w=0.2, h=0.3), alpha=0.9) +
  scale_color_manual(values = tmp.cols) +
  scale_size_manual(values = rep(1, length(tmp.cols))) +
  scale_shape_manual(values = rep(1, length(tmp.cols))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("coverage") +
  ylab("million CpG sites") +
  facet_wrap(~ downsampled_cov, ncol=n.cols.plot)


ggsave(filename= sprintf("CGsites_captured_%s.png", dat.name), 
       plot=p2, device= "png", width=10,height=h.plot, units="in") 
