# EPIQC CPG sites captured by Coverage plots
# Jessica Nordlund - June 2020

###------------###
#Load libs
###------------###
library(tidyr)
library(ggplot2)
library(reshape2)
options(stringsAsFactors = FALSE)

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
# Grab data: ie Cumulative coverage table
###------------###

setwd("/Users/jessicanordlund/Downloads/EPIQC")
test.data <- read.csv("cumulative_coverage_table_avgc_200414.txt", header= T, sep="\t") #Downsampled coverages

###---------------------###
#     SUBSET DATA 
###---------------------###

# 1- subset the data, biological replicates, remove genome merged and EPIC capture
dat <- test.data[!test.data$sample %in% "HG001-7" ,]
dat <- dat[!dat$prep %in% "EPIC" ,]

# 2 remove the 100 and 200x coverages and the 40x and 50x downsampling
dat <- dat[!dat$coverage_interval %in% c("[100,200)", "[200,Inf)") ,]
dat <- dat[!dat$downsampled %in% c("40x", "50x") ,]

# 3 Change the name of OXBS to TrueMethyl & then re-order alphabetically
dat$prep <- factor(dat$prep)
levels(dat$prep)[levels(dat$prep)=="OXBS"] <- "TrueMethyl"


# 4 Make numbers in Million w 2 digits
dat$value <- round(as.numeric(dat$cumulative_cg_site_count)/1000000, digits= 2)
dat$value2 <- round(as.numeric(dat$cg_site_count)/1000000, digits= 2)


# 4 Change lable names and order levels of factors
new.labels <- gsub("\\,.*", "", gsub("\\[", "", dat$coverage_interval))
dat$variable <- factor(new.labels, levels= unique(new.labels))
dat$downsampled_cov <- factor(dat$downsampled, levels = unique(dat$downsampled))
dat$variable2 <-  gsub(")", "]", dat$coverage_interval)
dat$variable2 <- factor(dat$variable2, levels = unique(dat$variable2))

# 5 Set plot colors by assay
tmp.cols <- my.assay.cols[match(levels(dat$prep), names(my.assay.cols))]

# Calc average number by assay
dat$variable3 <- paste(dat$prep, dat$downsampled, dat$variable2, sep="_")
unique.com <- unique(dat$variable3)

df <- data.frame()[1:length(unique.com), ]
for(i in unique.com){
  tmp.dat <- dat[dat$variable3 %in% i,]
  tmp.count <- c(i, 
                 unique(as.character(tmp.dat$prep)), 
                 unique(tmp.dat$downsampled),
                 unique(as.character(tmp.dat$variable)),
                 round(mean(tmp.dat$cumulative_cg_site_count), 0), 
                 round(mean(tmp.dat$value), 2))
  
  df <- rbind(df, tmp.count)
}   

names(df) <- c("Unique", "prep", "downsampled_cov", "variable", "CGs", "value")
df$value <- as.numeric(df$value)
df$variable <- factor(df$variable, levels=unique(df$variable))
df$downsampled_cov <- factor(df$downsampled_cov, levels=unique(df$downsampled_cov))

# Just getting some simple numbers for the manuscript text here
summary(dat[dat$downsampled %in% "10x" & dat$coverage_interval %in% "[5,10)" ,8])
summary(dat[dat$downsampled %in% "10x" & dat$coverage_interval %in% "[10,20)" ,8])
summary(dat[dat$downsampled %in% "10x" & dat$coverage_interval %in% "[5,10)" & dat$prep %in% "EMSeq",8])
summary(dat[dat$downsampled %in% "10x" & dat$coverage_interval %in% "[5,10)" & dat$prep %in% "TruSeq",8])

summary(dat[dat$downsampled %in% "20x" & dat$coverage_interval %in% "[10,20)" ,8])
summary(dat[dat$downsampled %in% "20x" & dat$coverage_interval %in% "[20,30)" ,8])

###---------------------###
#    Make the Plot
###---------------------###
p3 <- ggplot(df, aes(x = variable, y = value, fill = prep)) +
        geom_bar(stat = "identity", width=0.7, 
           position = 'dodge', 
           aes(color=prep), 
           alpha=0.7)+
      theme(aspect.ratio = 3/5)+
      scale_color_manual(values = tmp.cols) +
      scale_fill_manual(values = tmp.cols) +
      ylim(0,30) +
      theme_bw() +
      xlab("average coverage") +
      ylab("million CpG sites") +
      facet_wrap(~ downsampled_cov, ncol=2) +
      theme(legend.position="bottom") +
      theme(legend.title=element_blank()) +
      theme(text = element_text(size=12),
        axis.text.x = element_text(angle=0, size=12, hjust=0.5, vjust=0.5),
        panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))

ggsave(filename="AssayCompare_CGsites_captured.png", 
       plot=p3, device= "png", width=7,height=5, units="in") 

