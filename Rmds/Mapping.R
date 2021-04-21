# This Rmd generates the panels of Figure 1

# Dependencies
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(reshape2)
library(cowplot)
library(ggpubr)
library(data.table)
library(RColorBrewer)
library(ggExtra)
library(viridis)
library(ComplexHeatmap)

theme_linedraw2 = theme_linedraw() + theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
theme_jon = theme_linedraw2 + 
  theme(legend.position="bottom", 
        legend.title=element_text(size=11), 
        legend.text=element_text(size=10), 
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5), 
        axis.text.y=element_text(size=10), 
        axis.title=element_text(size=12), 
        axis.title.y=element_text(vjust=1), 
        plot.title = element_text(size=10, hjust = 0.5), 
        strip.text = element_text(size=12), 
        panel.grid.major = element_line(colour = "grey98"), 
        panel.grid.minor = element_blank())
ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}
my.assay.cols = c("MethylSeq" = "#377eb8", 
                  "TruSeq" = "#e41a1c",
                  "SPLAT" = "#984ea3",
                  "EMSeq" = "#4daf4a",
                  "EMSeqLAB01" = "#60c95d",
                  "EMSeqLAB02" = "#318f2e",
                  "EPIC" = "#d984c3",
                  "TrueMethyl" = "#a65628",
                  "TrueMethylBS" = "#a65628",
                  "TrueMethylOX" = "#fd8d3c",
                  "Nanopore" = "#4bcbde", "PromethION" = "#4bcbde")
my.assay.order = c("EMSeqLAB01", "EMSeqLAB02", "MethylSeq", "SPLAT", "TruSeq", "TrueMethylOX", "TrueMethylBS", "EPIC", "Nanopore")
my.assay2.order = c("EMSeq", "MethylSeq", "SPLAT", "TruSeq", "TrueMethylBS", "TrueMethylOX", "EPIC", "Nanopore")


# Mapping Rates (derived from samtools stats)
df_bwameth = read.csv("tables/samstats_all.csv")
df_bwameth = separate(df_bwameth, col="sampleOverall", sep="_", into=c("assay", "genome", "lab", "rep"), remove=F)
df_bwameth$multi2 = df_bwameth$multi + df_bwameth$nonprimary
df_bwameth$pctMapped = df_bwameth$mapped / df_bwameth$totalReads * 100
df_bwameth$pctMulti = df_bwameth$multi2 / df_bwameth$totalReads * 100
df_bwameth$pctDup = df_bwameth$duplicates / df_bwameth$totalReads * 100
df_bwameth$pctUnmapped = df_bwameth$unmapped / df_bwameth$totalReads * 100

# re-arrange names for coloring purposes
df_bwameth$assay2 = df_bwameth$assay
df_bwameth[which(df_bwameth$assay=="EMSeq" & df_bwameth$lab=="LAB01"), ]$assay = "EMSeqLAB01"
df_bwameth[which(df_bwameth$assay=="EMSeq" & df_bwameth$lab=="LAB02"), ]$assay = "EMSeqLAB02"
df_bwameth[which(df_bwameth$assay2=="TrueMethylOX"), ]$assay2 = "TrueMethyl"
df_bwameth[which(df_bwameth$assay2=="TrueMethylBS"), ]$assay2 = "TrueMethyl"
df_bwameth$assay = factor(df_bwameth$assay, levels=my.assay.order)
df_bwameth$assay2 = factor(df_bwameth$assay2, levels=my.assay2.order)

# Fig 1A: Mapping Efficiencies
gg_bwameth_map = ggplot(subset(df_bwameth, assay != "EPIC"), aes(x=assay2, y=pctMapped)) +
  geom_boxplot(aes(color=assay2, fill=assay2), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=assay), shape=21, show.legend=F, size=2) +
  scale_fill_manual(values=my.assay.cols) +
  scale_color_manual(values=my.assay.cols) +
  scale_x_discrete(labels = c("EMSeq"="EM", "MethylSeq"="MS", "SPLAT"="SP", "TruSeq"="TS", "TrueMethyl"="TM")) +
  theme_jon + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5)) +
  geom_hline(yintercept = 100) +
  xlab("") + ylab("% Read Pairs") +
  ggtitle("Primary Mapped")

gg_bwameth_multi = ggplot(subset(df_bwameth, assay != "EPIC"), aes(x=assay2, y=pctMulti)) +
  geom_boxplot(aes(color=assay2, fill=assay2), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=assay), shape=21, show.legend=F, size=2) +
  scale_fill_manual(values=my.assay.cols) +
  scale_color_manual(values=my.assay.cols) +
  scale_x_discrete(labels = c("EMSeq"="EM", "MethylSeq"="MS", "SPLAT"="SP", "TruSeq"="TS", "TrueMethyl"="TM")) +
  theme_jon + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5)) +
  xlab("") + ylab("") +
  ggtitle("Multi-Mapped")

gg_bwameth_dup = ggplot(subset(df_bwameth, assay != "EPIC"), aes(x=assay2, y=pctDup)) +
  geom_boxplot(aes(color=assay2, fill=assay2), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=assay), shape=21, show.legend=F, size=2) +
  scale_fill_manual(values=my.assay.cols) +
  scale_color_manual(values=my.assay.cols) +
  scale_x_discrete(labels = c("EMSeq"="EM", "MethylSeq"="MS", "SPLAT"="SP", "TruSeq"="TS", "TrueMethyl"="TM")) +
  theme_jon + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5)) +
  xlab("") + ylab("") +
  ggtitle("Duplicated")

gg_bwameth_unmapped = ggplot(subset(df_bwameth, assay != "EPIC"), aes(x=assay2, y=pctUnmapped)) +
  geom_boxplot(aes(color=assay2, fill=assay2), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=assay), shape=21, show.legend=F, size=2, cex = 0.3) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=my.assay.cols) +
  scale_color_manual(values=my.assay.cols) +
  scale_x_discrete(labels = c("EMSeq"="EM", "MethylSeq"="MS", "SPLAT"="SP", "TruSeq"="TS", "TrueMethyl"="TM")) +
  theme_jon + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5)) +
  xlab("") + ylab("") +
  ggtitle("Unmapped")

plot_grid(gg_bwameth_map, gg_bwameth_multi, gg_bwameth_dup, gg_bwameth_unmapped, rel_widths = c(3.2, 3, 3, 3), nrow=1)

# Bases Trimmed (derived from fastp reports)
df_fastp = read.csv("tables/fastp_report_all.csv")
df_fastp = df_fastp %>% group_by(sampleOverall) %>% summarise_at(.vars = colnames(df_fastp)[4:(ncol(df_fastp)-1)], .funs=sum)
df_fastp = separate(df_fastp, col="sampleOverall", sep="_", into=c("assay", "genome", "lab", "rep"), remove=F)
df_fastp$pctTrimmed = (df_fastp$prefilter_totalBases1 - df_fastp$postfilter_totalBases1) /  df_fastp$prefilter_totalBases1 * 100

df_fastp$assay2 = df_fastp$assay
df_fastp[which(df_fastp$assay2=="TrueMethylOX"), ]$assay2 = "TrueMethyl"
df_fastp[which(df_fastp$assay2=="TrueMethylBS"), ]$assay2 = "TrueMethyl"
df_fastp$assay = factor(df_fastp$assay, levels=c("EMSeq", "MethylSeq", "SPLAT", "TruSeq", "TrueMethylBS", "TrueMethylOX", "EPIC"))
df_fastp$assay2 = factor(df_fastp$assay2, levels=c("EMSeq", "MethylSeq", "SPLAT", "TruSeq", "TrueMethyl", "EPIC"))

gg_basestrimmed = ggplot(subset(df_fastp, assay != "EPIC"), aes(x=assay2, y=pctTrimmed)) +
  geom_boxplot(aes(color=assay2, fill=assay2), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=assay), shape=21, show.legend=F, size=2.5, cex=0.7) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(labels = c("EMSeq"="EM", "MethylSeq"="MS", "SPLAT"="SP", "TruSeq"="TS", "TrueMethyl"="TM")) +
  theme_jon + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5)) +
  scale_fill_manual(values=my.assay.cols) +
  scale_color_manual(values=my.assay.cols) +
  xlab("") + ylab("Percent") +
  ggtitle("Bases Trimmed")
gg_basestrimmed


# Insert Size (derived from samtools stats)
df_insert = read.csv("tables/samstats_insertsize.csv")
df_insert = separate(df_insert, col="sample", sep="_", into=c("assay", "genome", "lab", "rep"), remove=F)

# re-arrange names for coloring purposes
df_insert$assay2 = df_insert$assay
df_insert[which(df_insert$assay=="EMSeq" & df_insert$lab=="LAB01"), ]$assay = "EMSeqLAB01"
df_insert[which(df_insert$assay=="EMSeq" & df_insert$lab=="LAB02"), ]$assay = "EMSeqLAB02"
df_insert[which(df_insert$assay2=="TrueMethylOX"), ]$assay2 = "TrueMethyl"
df_insert[which(df_insert$assay2=="TrueMethylBS"), ]$assay2 = "TrueMethyl"
df_insert$assay = factor(df_insert$assay, levels=my.assay.order)
df_insert$assay2 = factor(df_insert$assay2, levels=c("EMSeq", "MethylSeq", "SPLAT", "TruSeq", "TrueMethyl", "EPIC"))

gg_insert = ggplot(subset(df_insert, assay != "EPIC"), aes(x=assay2, y=insertsize)) +
  geom_boxplot(aes(color=assay2, fill=assay2), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=assay), shape=21, show.legend=F, size=2.5, cex=0.7) +
  scale_fill_manual(values=my.assay.cols) +
  scale_color_manual(values=my.assay.cols) +
  scale_x_discrete(labels = c("EMSeq"="EM", "MethylSeq"="MS", "SPLAT"="SP", "TruSeq"="TS", "TrueMethyl"="TM")) +
  theme_jon + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5)) +
  xlab("") + ylab("Mean") +
  ggtitle("Insert Size")
gg_insert

# Cumulative Genome Coverage (derived from samtools stats)
df_cov = read.csv("tables/samstats_genomecov.csv")
df_cov = separate(df_cov, col="sample", sep="_", into=c("assay","genome"), remove=F)
df_cov[grepl("bs", df_cov$sample), ]$assay = "TrueMethylBS"
df_cov[grepl("ox", df_cov$sample), ]$assay = "TrueMethylOX"
df_cov = subset(df_cov, depth != 0)
df_cov = df_cov %>% group_by(sample) %>% mutate(cumulative_sum = cumsum(as.numeric(freq)))
df_cov$fraction_of_total = (df_cov$total - df_cov$cumulative_sum) / df_cov$total
df_cov

gg_genomecov = ggplot(df_cov, aes(x=depth, y=fraction_of_total*100)) +
  stat_smooth(size=2, span=0.1, se=F, aes(color=assay)) +
  geom_vline(xintercept = 20, linetype="dashed", alpha=0.6) +
  geom_hline(yintercept = 0, alpha=0.6) +
  theme_jon + theme(legend.position="none") +
  xlim(0, 200) +
  scale_color_manual(values = my.assay.cols) +
  xlab("Minimum Depth") + ylab("Percent") +
  ggtitle("Genomic Coverage")
gg_genomecov

# Dinucleotide Distribution (derived from Bismark reports)
df_dinuc = read.csv("tables/bismark_dinuc.csv")
df_dinuc = separate(df_dinuc, col="sample", sep="_", into=c("assay", "genome", "lab", "rep", "extra"), remove=F)
df_dinuc$difference = log2(df_dinuc$samplePct / df_dinuc$genomePct)
df_dinuc$dinucleotide = factor(df_dinuc$dinucleotide, levels=rev(c("C", "G", "CC", "CG", "GC", "GG", "CA", "CT",
                                                               "GA", "GT", "AC", "AG", "TC", "TG",
                                                               "A", "AA", "AT", "T", "TA", "TT")))

gg_dinuc = ggplot(df_dinuc, aes(x=dinucleotide, y=difference)) +
  geom_rect(xmin=0, xmax=6.5, ymin=-Inf, ymax=Inf, fill="#fff5f5", alpha=0.3) +
  geom_rect(xmin=14.5, xmax=Inf, ymin=-Inf, ymax=Inf, fill="#f0f6ff", alpha=0.3) +
  geom_line(aes(group=assay), alpha=0.5) +
  geom_point(aes(color=assay), size=2, show.legend=F) +
  geom_vline(linetype="dashed", xintercept = 6.5, alpha=0.3) +
  geom_vline(linetype="dashed", xintercept = 14.5, alpha=0.3) +
  geom_hline(linetype="dashed", yintercept = 0) +
  theme_jon + theme(legend.title = element_blank()) +
  scale_color_manual(values = my.assay.cols) +
  ylim(-0.75,0.75) +
  xlab("") + ylab("log2 Enrichment") +
  ggtitle("Dinucleotide Distribution") 
gg_dinuc


# CpG Depth versus Reads (and versus Bases) (derived from bwameth bedGraphs)
df_bwameth_CpGs = read.csv("tables/bwameth_CpGstats_all.csv")

# get sample name and mean depth
df_bwameth_CpGs_REP = df_bwameth_CpGs[, c("sample", "meanDepth")]
# get the total number of reads and bases
df_fastp = separate(df_fastp, col="sampleIndiv", sep="-", into=c("replevel", "rest"), remove=F)
df_fastp2 = df_fastp %>% group_by(replevel) %>% summarise(prefilter_totalReads1 = sum(prefilter_totalReads1), prefilter_totalBases1 = sum(prefilter_totalBases1))
colnames(df_fastp2) = c("sample", "totalReads", "totalBases")
df_bwameth_CpGs_REP = merge(df_bwameth_CpGs_REP, df_fastp2, by="sample")
df_bwameth_CpGs_REP = separate(df_bwameth_CpGs_REP, col="sample", sep="_", into=c("assay","genome","lab","rep"), remove=F)
# re-arrange names for coloring purposes
df_bwameth_CpGs_REP$assay2 = df_bwameth_CpGs_REP$assay
df_bwameth_CpGs_REP[which(df_bwameth_CpGs_REP$assay=="EMSeq" & df_bwameth_CpGs_REP$lab=="LAB01"), ]$assay = "EMSeqLAB01"
df_bwameth_CpGs_REP[which(df_bwameth_CpGs_REP$assay=="EMSeq" & df_bwameth_CpGs_REP$lab=="LAB02"), ]$assay = "EMSeqLAB02"
#df_bwameth_CpGs_REP[which(df_bwameth_CpGs_REP$assay2=="TrueMethylOX"), ]$assay2 = "TrueMethyl"
#df_bwameth_CpGs_REP[which(df_bwameth_CpGs_REP$assay2=="TrueMethylBS"), ]$assay2 = "TrueMethyl"

### manually add ONT
df_bwameth_CpGs_REP[nrow(df_bwameth_CpGs_REP)+1 , ] = c("Nanopore_HG001_LAB01_REP01", "Nanopore", "HG001", "LAB01", "REP01", 40.3905, NA,  77448216458, "Nanopore")
df_bwameth_CpGs_REP[nrow(df_bwameth_CpGs_REP)+1 , ] = c("Nanopore_HG002_LAB01_REP01", "Nanopore", "HG002", "LAB01", "REP01", 62.3843, NA, 224133506247, "Nanopore")
df_bwameth_CpGs_REP[nrow(df_bwameth_CpGs_REP)+1 , ] = c("Nanopore_HG003_LAB01_REP01", "Nanopore", "HG003", "LAB01", "REP01", 57.0771, NA, 210182000975, "Nanopore")
df_bwameth_CpGs_REP[nrow(df_bwameth_CpGs_REP)+1 , ] = c("Nanopore_HG004_LAB01_REP01", "Nanopore", "HG004", "LAB01", "REP01", 65.3639, NA, 232695284206, "Nanopore")
df_bwameth_CpGs_REP[nrow(df_bwameth_CpGs_REP)+1 , ] = c("Nanopore_HG005_LAB01_REP01", "Nanopore", "HG005", "LAB01", "REP01", 22.2183, NA, 113592270660, "Nanopore")
df_bwameth_CpGs_REP[nrow(df_bwameth_CpGs_REP)+1 , ] = c("Nanopore_HG006_LAB01_REP01", "Nanopore", "HG006", "LAB01", "REP01", 59.4287, NA, 209072335769, "Nanopore")
df_bwameth_CpGs_REP[nrow(df_bwameth_CpGs_REP)+1 , ] = c("Nanopore_HG007_LAB01_REP01", "Nanopore", "HG007", "LAB01", "REP01", 26.1968, NA, 102684455597, "Nanopore")
df_bwameth_CpGs_REP$meanDepth = as.double(df_bwameth_CpGs_REP$meanDepth)
df_bwameth_CpGs_REP$totalReads = as.double(df_bwameth_CpGs_REP$totalReads)
df_bwameth_CpGs_REP$totalBases = as.double(df_bwameth_CpGs_REP$totalBases)
df_bwameth_CpGs_REP$assay = factor(df_bwameth_CpGs_REP$assay, levels=my.assay.order)
df_bwameth_CpGs_REP$assay2 = factor(df_bwameth_CpGs_REP$assay2, levels=my.assay2.order)

gg_CpGs_vsReads = ggplot(subset(df_bwameth_CpGs_REP, assay != "EPIC" & sample != "TruSeq_HG006_LAB01_REP01" & sample != "TruSeq_HG007_LAB01_REP01" & assay2 != "Nanopore"), aes(y=meanDepth, x=totalReads)) +
  geom_hline(yintercept = 20, linetype="dashed", alpha=0.2) +
  geom_point(aes(fill=assay2), shape=21, position = "dodge", size=2.5, alpha=0.8) +
  geom_line(aes(color=assay2), stat="smooth", method = "lm", se=F, fullrange=TRUE, linetype="dashed", show.legend = F) +
  scale_color_manual(values=my.assay.cols) +
  scale_fill_manual(values=my.assay.cols) +
  theme_jon + theme(legend.position="none") +
  scale_x_continuous(breaks=c(100000000, 200000000, 300000000,400000000,500000000, 600000000, 700000000, 800000000), 
                     labels = c("100M", "200M", "300M", "400M","500M", "600M", "700M", "800M"), limits = c(50000000, 800000000)) +
  scale_y_continuous(breaks=c(0,10,20,30,40,50)) +
  ylab("Mean Depth per CpG") + xlab("Total Read Pairs") + 
  ggtitle("CpG Depth vs. Read Pairs") 
  
# plot with ONT
gg_CpGs_vsBases = ggplot(df_bwameth_CpGs_REP, aes(y=meanDepth, x=totalBases)) +
  geom_point(data=subset(df_bwameth_CpGs_REP, assay2 != "Nanopore"), aes(fill=assay2), shape=21, position = "dodge", size=2.5, alpha=0.15) +
  geom_point(data=subset(df_bwameth_CpGs_REP, assay2 == "Nanopore"), aes(fill=assay2), shape=21, position = "dodge", size=2.5) +
  geom_line(data=subset(df_bwameth_CpGs_REP, assay2 != "Nanopore"), aes(color=assay2), stat="smooth", method = "lm", se=F, fullrange=TRUE, linetype="dashed", alpha=0.15, show.legend = F) +
  geom_line(data=subset(df_bwameth_CpGs_REP, assay2 == "Nanopore"), aes(color=assay2), stat="smooth", method = "lm", se=F, fullrange=TRUE, linetype="dashed", show.legend = F) +
  scale_color_manual(values=my.assay.cols) +
  scale_fill_manual(values=my.assay.cols) +
  theme_jon + theme(legend.position="right", legend.title=element_blank()) +
  scale_x_continuous(breaks=c(50000000000, 100000000000, 150000000000, 200000000000), 
                     labels = c("500M", "1.0B", "1.5B", "2.0B")) +
  scale_y_continuous(breaks=seq(0,100,10)) +
  ylab("") + xlab("Total Bases") + 
  ggtitle("CpG Depth vs. Total Bases") 
gg_CpGs_vsBases
