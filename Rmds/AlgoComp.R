EpiQC // Algorithm Comparison

# Dependencies
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(reshape2)
library(cowplot)
library(ggrepel)
library(data.table)
library(viridis)
library(ggpubr)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

theme_linedraw2 = theme_linedraw() + theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
theme_jon = theme_linedraw2 + 
  theme(legend.position="bottom", 
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12), 
        axis.text.x = element_text(size=12, angle=90, hjust=1, vjust=0.5), 
        axis.text.y=element_text(size=12), 
        axis.title=element_text(size=14), 
        axis.title.y=element_text(vjust=1, size=11),
        axis.title.x = element_text(size=11),
        plot.title = element_text(size=13, hjust = 0.5), 
        strip.text = element_text(size=12), 
        panel.grid.major = element_line(colour = "grey98"), 
        panel.grid.minor = element_blank())
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
my.assay.order = c("EMSeqLAB01", "EMSeqLAB02", "MethylSeq", "SPLAT", "TruSeq", "TrueMethylOX", "TrueMethylBS", "EPIC")
my.assay2.order = c("EMSeq", "MethylSeq", "SPLAT", "TruSeq", "TrueMethyl", "EPIC")
my.algo.cols = c("bismark" = "#fdc086", "bwameth"="#386cb0", "bsseeker2"="#7fc97f", "bitmapperbs"="#ed6baf", "bitmapperBS"="#ed6baf")
my.algo.cols2 = c("bis" = "#fdc086", "bwa"="#386cb0", "bsk"="#7fc97f", "bit"="#ed6baf")

# Cross-Algorithm Mapping Comparison
df_algo_map = read.csv("tables/samstats_allalgos.csv")
df_algo_map = separate(df_algo_map, col="sampleOverall", sep="_", into=c("assay", "genome", "lab", "rep"), remove=F)
df_algo_map$pctMapped = ifelse(df_algo_map$algo == "bitmapperBS", (df_algo_map$mapped - df_algo_map$unmapped) / df_algo_map$totalReads * 100 , df_algo_map$mapped / df_algo_map$totalReads * 100)
df_algo_map$pctMulti = df_algo_map$multi / df_algo_map$totalReads * 100
df_algo_map$pctDup = df_algo_map$duplicates / df_algo_map$totalReads * 100
df_algo_map$pctUnmapped = df_algo_map$unmapped / df_algo_map$totalReads * 100

# re-arrange names for coloring purposes
df_algo_map$assay2 = df_algo_map$assay
df_algo_map[which(df_algo_map$assay=="EMSeq" & df_algo_map$lab=="LAB01"), ]$assay = "EMSeqLAB01"
df_algo_map[which(df_algo_map$assay=="EMSeq" & df_algo_map$lab=="LAB02"), ]$assay = "EMSeqLAB02"
df_algo_map[which(df_algo_map$assay2=="TrueMethylOX"), ]$assay2 = "TrueMethyl"
df_algo_map[which(df_algo_map$assay2=="TrueMethylBS"), ]$assay2 = "TrueMethyl"
df_algo_map$assay = factor(df_algo_map$assay, levels=my.assay.order)
df_algo_map$assay2 = factor(df_algo_map$assay2, levels=my.assay2.order)

# Mapping Efficiency Plot
gg_algo_map = ggplot(df_algo_map, aes(x=algo, y=pctMapped)) +
  geom_boxplot(aes(color=algo, fill=algo), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=algo), shape=21, show.legend=F, size=2) +
  scale_fill_manual(values=my.algo.cols) +
  scale_color_manual(values=my.algo.cols) +
  scale_x_discrete(labels = c("bismark"="Bis", "bwameth"="Bwa", "bsseeker2"="Bsk", "bitmapperBS"="Bit")) +
  theme_jon + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5, size=10), plot.title = element_text(size=11)) +
  geom_signif() +
  geom_hline(yintercept = 100) +
  xlab("") + ylab("% Read Pairs") +
  ggtitle("Primary Mapped")

gg_algo_multi = ggplot(subset(df_algo_map, pctMulti < 20), aes(x=algo, y=pctMulti)) +
  geom_boxplot(aes(color=algo, fill=algo), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=algo), shape=21, show.legend=F, size=2) +
  scale_fill_manual(values=my.algo.cols) +
  scale_color_manual(values=my.algo.cols) +
  scale_x_discrete(labels = c("bismark"="Bis", "bwameth"="Bwa", "bsseeker2"="Bsk", "bitmapperBS"="Bit")) +
  theme_jon + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5, size=10), plot.title = element_text(size=11)) +
  geom_hline(yintercept = 0) +
  xlab("") + ylab("") +
  ggtitle("Multi-Mapped")

gg_algo_dup = ggplot(subset(df_algo_map, assay != "EPIC"), aes(x=algo, y=pctDup)) +
  geom_boxplot(aes(color=algo, fill=algo), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=algo), shape=21, show.legend=F, size=2) +
  scale_fill_manual(values=my.algo.cols) +
  scale_color_manual(values=my.algo.cols) +
  scale_x_discrete(labels = c("bismark"="Bis", "bwameth"="Bwa", "bsseeker2"="Bsk", "bitmapperBS"="Bit")) +
  theme_jon + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5, size=10), plot.title = element_text(size=11)) +
  xlab("") + ylab("") +
  ggtitle("Duplicated")

gg_algo_unmapped = ggplot(df_algo_map, aes(x=algo, y=pctUnmapped)) +
  geom_boxplot(aes(color=algo, fill=algo), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=algo), shape=21, show.legend=F, size=2, cex = 0.3) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=my.algo.cols) +
  scale_color_manual(values=my.algo.cols) +
  scale_x_discrete(labels = c("bismark"="Bis", "bwameth"="Bwa", "bsseeker2"="Bsk", "bitmapperBS"="Bit")) +
  theme_jon + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5, size=10), plot.title = element_text(size=11)) +
  xlab("") + ylab("") +
  ggtitle("Unmapped")

gg_algo_all = plot_grid( gg_algo_map, gg_algo_multi, gg_algo_dup, gg_algo_unmapped, rel_widths = c(3.3, 3, 3, 3.1), nrow=1)
gg_algo_all

# Bases Mapped and CpG Depth
df_algo_efficiencyAndCpG = read.csv("tables/allalgos_basesMappedAndCpGdepth.csv")
df_algo_efficiencyAndCpG = separate(df_algo_efficiencyAndCpG, col="sampleOverall", into=c("assay", "genome", "lab", "rep"), remove=F)
df_algo_efficiencyAndCpG$efficiency = df_algo_efficiencyAndCpG$mapped_bases / df_algo_efficiencyAndCpG$total_bases * 100
df_algo_efficiencyAndCpG[df_algo_efficiencyAndCpG$assay=="TrueMethylOX", ]$assay = "TrueMethyl"
df_algo_efficiencyAndCpG[df_algo_efficiencyAndCpG$assay=="TrueMethylBS", ]$assay = "TrueMethyl"
df_algo_efficiencyAndCpG = df_algo_efficiencyAndCpG %>% group_by(sampleOverall, algo, assay) %>% summarise(CpGdepth = sum(CpGdepth), efficiency = mean(efficiency))
df_algo_efficiencyAndCpG = separate(df_algo_efficiencyAndCpG, col="sampleOverall", sep="_", into=c("assay","genome","lab","rep"), remove=F)
df_algo_efficiencyAndCpG$labrep = paste0(df_algo_efficiencyAndCpG$lab, "\n", df_algo_efficiencyAndCpG$rep)
df_algo_efficiencyAndCpG[df_algo_efficiencyAndCpG$assay=="TrueMethylOX", ]$assay = "TrueMethyl"
df_algo_efficiencyAndCpG[df_algo_efficiencyAndCpG$assay=="TrueMethylBS", ]$assay = "TrueMethyl"
df_algo_efficiencyAndCpG$algo = factor(df_algo_efficiencyAndCpG$algo, levels=c("bismark", "bsseeker2", "bitmapperbs", "bwameth"))

sample.labels = c("EMSeq_HG002_LAB01_REP01" = "L1R1",
"EMSeq_HG002_LAB01_REP02" = "L1R2",
"EMSeq_HG002_LAB02_REP01" = "L2R1",
"EMSeq_HG002_LAB02_REP02" = "L2R2",
"MethylSeq_HG002_LAB01_REP01" = "L1R1",
"MethylSeq_HG002_LAB01_REP02" = "L1R2",
"SPLAT_HG002_LAB01_REP01" = "L1R1",
"SPLAT_HG002_LAB01_REP02" = "L1R2",
"TrueMethylBS_HG002_LAB01_REP01" = "BS1",
"TrueMethylBS_HG002_LAB01_REP02" = "BS2",
"TrueMethylOX_HG002_LAB01_REP01" = "OX1",
"TrueMethylOX_HG002_LAB01_REP02" = "OX2",
"TruSeq_HG002_LAB01_REP01" = "L1R1",
"TruSeq_HG002_LAB01_REP02" = "L1R2")

gg_efficiency = ggplot(df_algo_efficiencyAndCpG, aes(x=sampleOverall, y=efficiency, group=algo)) +
  geom_hline(yintercept = 100, alpha=0.5) +
  geom_linerange(aes(x=sampleOverall, ymin=60, ymax=efficiency), alpha=0.5, position = position_dodge(width=0.6)) +
  geom_point(aes(fill=algo), shape=21, show.legend = F, size=2.75, position = position_dodge(width=0.6)) +
  scale_fill_manual(values=my.algo.cols) +
  scale_x_discrete(labels = sample.labels) +
  theme_jon + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5, size=10)) +
  xlab("") + ylab("Percent") + ggtitle("Bases Mapped") +
  facet_wrap(~assay, nrow=1, scales="free_x")
gg_efficiency


df_algo_efficiencyAndCpG = df_algo_efficiencyAndCpG %>% group_by(sampleOverall, algo, assay) %>% summarise(CpGdepth = sum(CpGdepth), efficiency = mean(efficiency))
df_algo_efficiencyAndCpG = separate(df_algo_efficiencyAndCpG, col="sampleOverall", sep="_", into=c("assay","genome","lab","rep"), remove=F)
df_algo_efficiencyAndCpG$labrep = paste0(df_algo_efficiencyAndCpG$lab, "\n", df_algo_efficiencyAndCpG$rep)
df_algo_efficiencyAndCpG[df_algo_efficiencyAndCpG$assay=="TrueMethylOX", ]$assay = "TrueMethyl"
df_algo_efficiencyAndCpG[df_algo_efficiencyAndCpG$assay=="TrueMethylBS", ]$assay = "TrueMethyl"
df_algo_efficiencyAndCpG$algo = factor(df_algo_efficiencyAndCpG$algo, levels=c("bismark", "bsseeker2", "bitmapperbs", "bwameth"))

gg_CpGdepth = ggplot(df_algo_efficiencyAndCpG, aes(x=sampleOverall, y=CpGdepth, group=algo)) +
  geom_hline(yintercept = 0, alpha=0.5) +
  geom_linerange(aes(x=sampleOverall, ymin=0, ymax=CpGdepth), alpha=0.5, position = position_dodge(width=0.6)) +
  geom_point(aes(fill=algo), shape=21, show.legend = T, size=2.75, position = position_dodge(width=0.6)) +
  scale_fill_manual(values=my.algo.cols) +
  scale_x_discrete(labels = sample.labels) +
  scale_y_continuous(breaks = c(0,5,10,15)) + coord_cartesian(ylim=c(0,15)) +
  theme_jon + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5, size=10), legend.title = element_blank()) +
  xlab("") + ylab("Mean") + ggtitle("CpG Coverage") +
  facet_wrap(~assay, nrow=1, scales="free_x")
gg_CpGdepth



# Annotate genome capture of ds20
results = data.frame()
bgs = list.files(path="bedgraphs/", pattern="chr1.ds20")
for(bg in bgs){
  curbg = read.csv(paste0("bedgraphs/", bg), sep="\t", header = F)  
  colnames(curbg) = c("chr", "start", "end", "meth", "numC", "numT")  
  
  # subsample for testing
  # curbg = sample_n(curbg, 100000)
    
  # get e.g. "EMSeq.bwa"
  compname = paste0(strsplit(bg, '_')[[1]][1], ".", strsplit(bg, '\\.')[[1]][4])
  curbg$Comparison = compname
  curbg$locus = paste0(curbg$chr, '-', curbg$start)
  results = rbind(results, curbg)
}
GR_results = makeGRangesFromDataFrame(results, seqnames.field = "chr", start.field = "start", end.field = "end",  keep.extra.columns = T)

geneRegionAnn = genomation::readTranscriptFeatures("hg38.refseqGenes.bed")

annots = c("hg38_cpgs", "hg38_enhancers_fantom", "hg38_genes_1to5kb", "hg38_genes_promoters", 
           "hg38_genes_5UTRs", "hg38_genes_cds", "hg38_genes_introns", "hg38_genes_3UTRs")
annotations = build_annotations(genome = 'hg38', annotations = annots)
peaksAnnotated = annotate_regions(
	regions = GR_results,
	annotations = annotations,
	ignore.strand = TRUE,
	quiet = FALSE,
	minoverlap = 1)
peaksAnnotated = data.frame(peaksAnnotated)
peaksAnnotated$locus = paste0(peaksAnnotated$seqnames, '-', peaksAnnotated$start)

## summarize each into different categories
## using precedence assign only one subtype per peak, per annotation category
peaksAnnotatedGene = peaksAnnotated[ grep("gene", peaksAnnotated$annot.type), ]
peaksAnnotatedGene$GenicAnnotation = mapvalues(peaksAnnotatedGene$annot.type, 
													c("hg38_genes_1to5kb", "hg38_genes_introns", "hg38_genes_3UTRs", "hg38_genes_5UTRs", "hg38_genes_promoters", "hg38_genes_cds"), 
													c("Up5kb", "Intron", "3'UTR", "5'UTR", "Promoter1kb", "CDS"))
peaksAnnotatedGene$GenicAnnotation = factor(peaksAnnotatedGene$GenicAnnotation, levels = c("Up5kb", "Intron", "Promoter1kb", "3'UTR", "5'UTR", "CDS"))
peaksAnnotatedGene = peaksAnnotatedGene[ order(peaksAnnotatedGene$GenicAnnotation, decreasing = T), ]
peaksAnnotatedGene = peaksAnnotatedGene[ !duplicated(peaksAnnotatedGene$locus), ]
peaksAnnotatedGene = peaksAnnotatedGene[, c("locus", "annot.symbol", "GenicAnnotation")]
colnames(peaksAnnotatedGene) = c("locus", "AnnotationGeneSymbol", "GenicAnnotation")
peaksAnnotatedGene$GenicAnnotation = factor(peaksAnnotatedGene$GenicAnnotation, levels = c("Intergenic", "Up5kb", "Promoter1kb", "5'UTR", "CDS", "Intron", "3'UTR"))
peaksAnnotatedCpG = peaksAnnotated[ grep("cpg", peaksAnnotated$annot.type), ]
peaksAnnotatedCpG$CpGAnnotation = mapvalues(peaksAnnotatedCpG$annot.type, 
													c("hg38_cpg_inter", "hg38_cpg_shelves", "hg38_cpg_shores", "hg38_cpg_islands"), 
													c("Intergenic", "CpG shelves", "CpG shores", "CpG islands"))
peaksAnnotatedCpG$CpGAnnotation = factor(peaksAnnotatedCpG$CpGAnnotation, levels = c("Intergenic", "CpG shelves", "CpG shores", "CpG islands"))
peaksAnnotatedCpG = peaksAnnotatedCpG[ order(peaksAnnotatedCpG$CpGAnnotation, decreasing = T), ]
peaksAnnotatedCpG = peaksAnnotatedCpG[ !duplicated(peaksAnnotatedCpG$locus), ]
peaksAnnotatedCpG = peaksAnnotatedCpG[, c("locus", "CpGAnnotation")]
peaksAnnotatedCpG$CpGAnnotation = factor(peaksAnnotatedCpG$CpGAnnotation, levels = c("CpG islands", "CpG shores", "CpG shelves", "Intergenic"))
peaksAnnotatedFantom = peaksAnnotated[ grep("enhancers", peaksAnnotated$annot.type), ]
peaksAnnotatedFantom$FantomEnhancerAnnotation = mapvalues(peaksAnnotatedFantom$annot.type, 
													c("hg38_enhancers_fantom"), 
													c("Enhancer"))
peaksAnnotatedFantom$FantomEnhancerAnnotation = factor(peaksAnnotatedFantom$FantomEnhancerAnnotation, levels = c("Enhancer"))
peaksAnnotatedFantom = peaksAnnotatedFantom[ order(peaksAnnotatedFantom$FantomEnhancerAnnotation, decreasing = T), ]
peaksAnnotatedFantom = peaksAnnotatedFantom[ !duplicated(peaksAnnotatedFantom$locus), ]
peaksAnnotatedFantom = peaksAnnotatedFantom[, c("locus", "FantomEnhancerAnnotation")]
peaksAnnotatedFantom$FantomEnhancerAnnotation = factor(peaksAnnotatedFantom$FantomEnhancerAnnotation, levels = c("Enhancer", "Non-enhancer"))
peaksAnnotatedSum = results[, c("locus", "chr", "start", "end")]
peaksAnnotatedSum = merge(peaksAnnotatedSum, peaksAnnotatedGene, by="locus", all.x = T)
peaksAnnotatedSum = merge(peaksAnnotatedSum, peaksAnnotatedCpG, by="locus", all.x = T)
peaksAnnotatedSum = merge(peaksAnnotatedSum, peaksAnnotatedFantom, by="locus", all.x = T)
peaksAnnotatedSum$GenicAnnotation[ is.na(peaksAnnotatedSum$GenicAnnotation) ] = "Intergenic"
peaksAnnotatedSum$FantomEnhancerAnnotation[ is.na(peaksAnnotatedSum$FantomEnhancerAnnotation) ] = "Non-enhancer"

resultsAnn = merge(results, peaksAnnotatedSum[, -c(2:4)], by = "locus") 

resultsSum = data.frame()
for(curComp in unique(results$Comparison))
{
	curRes = resultsAnn[ resultsAnn$Comparison %in% curComp, ]
	curSum = curRes %>% group_by(GenicAnnotation) %>% dplyr::summarise(Count = n())
	colnames(curSum)[1] = "Annotation"
	resultsSum = rbind(resultsSum, data.frame(Comparison = curComp, AnnotationType = "Genic", curSum))
	
	curSum = curRes %>% group_by(CpGAnnotation) %>% dplyr::summarise(Count = n())
	colnames(curSum)[1] = "Annotation"
	resultsSum = rbind(resultsSum, data.frame(Comparison = curComp, AnnotationType = "CpG", curSum))
	
	curSum = curRes %>% group_by(FantomEnhancerAnnotation) %>% dplyr::summarise(Count = n())
	colnames(curSum)[1] = "Annotation"
	resultsSum = rbind(resultsSum, data.frame(Comparison = curComp, AnnotationType = "Fantom Enhancers", curSum))
}
resultsSum2 = resultsSum

annTypeColors = c("Intergenic"="#DDDEE0", "Up5kb"="#F49014", "Promoter1kb"="#DF6A0B", "5'UTR"="#134F86", "CDS"="#06882C", "Intron"="#fbb4ae", "3'UTR"="#763F9A", 
  "CpG shelves"="#55AA55", "CpG shores"="#2D882D", "CpG islands"="#004400", 
  "Non-enhancer"="#DDDEE0", "Enhancer"="#e7298a")
resultsSum2$Annotation = factor(resultsSum2$Annotation, levels = names(annTypeColors))
resultsSum2 = separate(resultsSum2, col="Comparison", sep="\\.", into=c("assay", "protocol"), remove=F)

ggplot(subset(resultsSum2, AnnotationType == "Genic"), 
       aes(x=assay, y=Count, fill=protocol)) + 
  geom_linerange(aes(x=assay, ymin=0, ymax=Count), alpha=0.5, position = position_dodge(width=0.6)) +
  geom_point(aes(fill=protocol), shape=21, show.legend = F, size=2.75, position = position_dodge(width=0.6)) +
  theme_jon +
  xlab("") + ylab("Count") +
  facet_wrap(~Annotation, scales = "free_y", nrow=1) +
  scale_fill_manual(values=my.algo.cols2)
