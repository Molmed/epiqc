# Dependencies
library(dplyr)
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
library(GenomicRanges)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
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

# CpG Coverage on Chromosome 1
rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
r  = rf(32)

df_cpgcov = read.csv("bedgraphs/cormatrix/corMatrix_coverage_ds20.chr1.csv.gz")
df_cpgcov0 = df_cpgcov#[1:500, ] for testing
df_cpgcov0[2:ncol(df_cpgcov0)] = lapply(df_cpgcov0[2:ncol(df_cpgcov0)], function(x) ifelse(is.na(x), 0, ifelse(x < 5, 0, 1)))
df_cpgcov2 = df_cpgcov0[, 2:ncol(df_cpgcov0)]
df_cpgcov2 = df_cpgcov2[do.call(order, as.list(df_cpgcov2)), ]
df_cpgcov2$id = 1:nrow(df_cpgcov2)
df_cpgcov3 = melt(df_cpgcov2, "id")
df_cpgcov3$value = factor(df_cpgcov3$value, levels=c(1, 0))

gg_altuna = ggplot(df_cpgcov3, aes(x=rev(id), y=variable, fill=value, color=value)) + 
          geom_tile(show.legend=T) + scale_color_manual(values=c("#fed976", "#0570b0"), labels=c("Present", "<5x Cov")) + 
          scale_fill_manual(values=c("#fed976", "#0570b0"), labels=c("Present", "<5x Cov")) +  
          xlab("CpGs") + ylab("") + theme_jon + theme(legend.title = element_blank()) + ggtitle("Chr1 CpG Coverage (20X Downsample)")
gg_altuna

# CpG annotation, cov, & meth
results = data.frame()
bgs = list.files(path="bedgraphs/", pattern="chr1.ds20.bwa")
for(bg in bgs){
  curbg = read.csv(paste0("bedgraphs/", bg), sep="\t", header = F)  
  colnames(curbg) = c("chr", "start", "end", "meth", "numC", "numT")  
  
  # subsample for testing
  # curbg = sample_n(curbg, 100000)
  
  # filter for min 10x
  curbg$depth = curbg$numC + curbg$numT
  #curbg = subset(curbg, depth > 19)
  
  # get assay name
  compname = strsplit(bg, '_')[[1]][1]
  curbg$Comparison = compname
  curbg$locus = paste0(curbg$chr, '-', curbg$start)
  results = rbind(results, curbg)
}


GR_results = makeGRangesFromDataFrame(results, seqnames.field = "chr", start.field = "start", end.field = "end",  keep.extra.columns = T)
geneRegionAnn = genomation::readTranscriptFeatures("hg38.refseqGenes.bed")
GR_cpgIslands = getGR_CpGi("hg38.cpgi.bed")
CpGisland = data.frame()
for(comparison in unique(results$Comparison)){
  cur_CpGisland = annotateCpGi(subset(GR_results, Comparison == comparison), GR_cpgIslands)
  cur_CpGisland$Comparison = comparison
  CpGisland = rbind(CpGisland, cur_CpGisland)
}
CpGisland$CountTotal = ave(CpGisland$Count, CpGisland$Comparison, FUN = sum)
CpGisland$CountPct = CpGisland$Count / CpGisland$CountTotal * 100

codingRegionAnn_all = data.frame()
for(comparison in unique(results$Comparison)){
  cur_codingRegionAnn = data.frame(genomation::getTargetAnnotationStats(genomation::annotateWithGeneParts(subset(GR_results, Comparison == comparison), geneRegionAnn), percentage=FALSE, precedence=TRUE))
  cur_codingRegionAnn$Comparison = comparison
  codingRegionAnn_all = rbind(codingRegionAnn_all, cur_codingRegionAnn)
}
colnames(codingRegionAnn_all) = c("Count", "Comparison")
codingRegionAnn_all$region = rownames(codingRegionAnn_all)
codingRegionAnn_all$region = gsub("[[:digit:]]+", '', codingRegionAnn_all$region)

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
resultsAnn$Comparison = gsub("PromethION", "Nanopore", resultsAnn$Comparison)

resultsSum = data.frame()
for(curComp in unique(results$Comparison))
{
	if(curComp == "PromethION") curComp = "Nanopore"
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
resultsSum2

annTypeColors = c("Intergenic"="#DDDEE0", "Up5kb"="#F49014", "Promoter1kb"="#DF6A0B", "5'UTR"="#134F86", "CDS"="#06882C", "Intron"="#fbb4ae", "3'UTR"="#763F9A", 
  "CpG shelves"="#55AA55", "CpG shores"="#2D882D", "CpG islands"="#004400", 
  "Non-enhancer"="#DDDEE0", "Enhancer"="#e7298a")
resultsSum2$Annotation = factor(resultsSum2$Annotation, levels = rev(c("Up5kb", "Promoter1kb", "5'UTR", "CDS", "Intron", "3'UTR", "Intergenic", "CpG shelves", "CpG shores", "CpG islands", "Enhancer", "Non-enhancer")))


# Global genomic coverage
gg_CpGcovGenomewide = ggplot(subset(resultsSum2, AnnotationType=="Genic"), 
       aes(x=Comparison, y=Count, fill=Annotation)) + 
  #geom_linerange(aes(x=assay, ymin=0, ymax=Count), alpha=0.5, position = position_dodge(width=0.6)) +
  #geom_point(aes(fill=protocol), shape=21, show.legend = F, size=2.75, position = position_dodge(width=0.6)) +
  geom_bar(stat="identity", position="fill") +
  theme_jon + theme(legend.title = element_blank(), legend.position="right", legend.text=element_text(size=10)) +
  xlab("") + ylab("Fraction") + ggtitle("All CpGs") + coord_flip() +
  scale_fill_manual(values=annTypeColors, labels=rev(c("Up5kb", "Promoter", "5'UTR", "Exon", "Intron", "3'UTR", "Intergenic"))) 


## Get distributions of UNIQUELY CAUGHT (only one platform got this CpG)
resultsSumSub = data.frame()
for(plat in unique(results$Comparison)) {
  if(plat=="PromethION") plat = "Nanopore"
  curPlatUniqLoci = df_cpgcov0[ df_cpgcov0[[plat]]==1 & rowSums(df_cpgcov0[,2:7])==1, ]$locus
  # stupid hack
  curPlatUniqLoci = gsub("chr1.","", curPlatUniqLoci)
  curPlatUniqLoci = gsub("\\..*", "", curPlatUniqLoci)
  curPlatUniqLoci = paste0("chr1-", curPlatUniqLoci)
  
  # subset resultsAnn
  resultsAnnSub = subset(resultsAnn, locus %in% curPlatUniqLoci & Comparison == plat)
  if(nrow(resultsAnnSub) == 0) resultsAnnSub = head(subset(resultsAnn, Comparison == plat), 1)
  
  curComp = plat
  curRes = resultsAnnSub[ resultsAnnSub$Comparison %in% curComp, ]
  curSum = curRes %>% group_by(GenicAnnotation) %>% dplyr::summarise(Count = n())
  colnames(curSum)[1] = "Annotation"
  resultsSumSub = rbind(resultsSumSub, data.frame(Comparison = curComp, AnnotationType = "Genic", curSum))
  
  curSum = curRes %>% group_by(CpGAnnotation) %>% dplyr::summarise(Count = n())
  colnames(curSum)[1] = "Annotation"
  resultsSumSub = rbind(resultsSumSub, data.frame(Comparison = curComp, AnnotationType = "CpG", curSum))
  
  curSum = curRes %>% group_by(FantomEnhancerAnnotation) %>% dplyr::summarise(Count = n())
  colnames(curSum)[1] = "Annotation"
  resultsSumSub = rbind(resultsSumSub, data.frame(Comparison = curComp, AnnotationType = "Fantom Enhancers", curSum))
}

resultsSumSub$Annotation = factor(resultsSumSub$Annotation, levels = rev(c("Up5kb", "Promoter1kb", "5'UTR", "CDS", "Intron", "3'UTR",
                                                                               "Intergenic", "CpG shelves", "CpG shores", "CpG islands",
                                                                               "Enhancer", "Non-enhancer")))

title = ggdraw() + draw_label("Uniquely Covered", x=0.3, hjust=0) + theme(plot.margin=margin(0,0,0,0))
gg_uniquelycovered = plot_grid(
  
  title,
  
  plot_grid(ggplot(subset(resultsSumSub, AnnotationType=="Genic"), 
       aes(x=Comparison, y=Count+50)) + 
  geom_bar(stat="identity", position="stack") +
  theme_jon + theme(legend.title = element_blank()) +
  xlab("") + ylab("Count") +  coord_flip() +
  scale_fill_manual(values=annTypeColors),
  
  ggplot(subset(resultsSumSub, AnnotationType=="Genic"), 
       aes(x=Comparison, y=Count, fill=Annotation)) + 
  geom_bar(stat="identity", position="fill") +
  theme_jon + theme(legend.title = element_blank(), legend.position="none", axis.text.y = element_blank()) +
  xlab("") + ylab("Fraction") + coord_flip() +
  scale_fill_manual(values=annTypeColors), nrow=1, rel_widths=c(0.5,0.7)),
  nrow=2, rel_heights=c(0.1,0.9)
)
gg_uniquelycovered


## Get distributions of UNIQUELY MISSED (only one platform didn't get this CpG)
resultsSumSub = data.frame()
for(plat in unique(results$Comparison)) {
  if(plat=="PromethION") plat = "Nanopore"
  curPlatUniqLoci = df_cpgcov0[ df_cpgcov0[[plat]]==0 & rowSums(df_cpgcov0[,2:7])==5, ]$locus
  # stupid hack
  curPlatUniqLoci = gsub("chr1.","", curPlatUniqLoci)
  curPlatUniqLoci = gsub("\\..*", "", curPlatUniqLoci)
  curPlatUniqLoci = paste0("chr1-", curPlatUniqLoci)
  
  # subset resultsAnn
  resultsAnnSub = subset(resultsAnn, locus %in% curPlatUniqLoci & Comparison == plat)
  if(nrow(resultsAnnSub) == 0) resultsAnnSub = head(subset(resultsAnn, Comparison == plat), 1)
  
  curComp = plat
  curRes = resultsAnnSub[ resultsAnnSub$Comparison %in% curComp, ]
  curSum = curRes %>% group_by(GenicAnnotation) %>% dplyr::summarise(Count = n())
  colnames(curSum)[1] = "Annotation"
  resultsSumSub = rbind(resultsSumSub, data.frame(Comparison = curComp, AnnotationType = "Genic", curSum))
  
  curSum = curRes %>% group_by(CpGAnnotation) %>% dplyr::summarise(Count = n())
  colnames(curSum)[1] = "Annotation"
  resultsSumSub = rbind(resultsSumSub, data.frame(Comparison = curComp, AnnotationType = "CpG", curSum))
  
  curSum = curRes %>% group_by(FantomEnhancerAnnotation) %>% dplyr::summarise(Count = n())
  colnames(curSum)[1] = "Annotation"
  resultsSumSub = rbind(resultsSumSub, data.frame(Comparison = curComp, AnnotationType = "Fantom Enhancers", curSum))
}

resultsSumSub$Annotation = factor(resultsSumSub$Annotation, levels = rev(c("Up5kb", "Promoter1kb", "5'UTR", "CDS", "Intron", "3'UTR",
                                                                               "Intergenic", "CpG shelves", "CpG shores", "CpG islands",
                                                                               "Enhancer", "Non-enhancer")))

title = ggdraw() + draw_label("Uniquely Missed", x=0.15, hjust=0) + theme(plot.margin=margin(0,0,0,0))
gg_uniquelymissed = plot_grid(
  
  title,
  
  plot_grid(ggplot(subset(resultsSumSub, AnnotationType=="Genic"), 
       aes(x=Comparison, y=Count+50)) + 
  geom_bar(stat="identity", position="stack") +
  theme_jon + theme(legend.title = element_blank(), axis.text.y = element_blank()) +
  xlab("") + ylab("Count") +  coord_flip() +
  scale_fill_manual(values=annTypeColors),
  
  ggplot(subset(resultsSumSub, AnnotationType=="Genic"), 
       aes(x=Comparison, y=Count, fill=Annotation)) + 
  geom_bar(stat="identity", position="fill") +
  theme_jon + theme(legend.title = element_blank(), legend.position="right", legend.text = element_text(size=9), axis.text.y = element_blank()) +
  xlab("") + ylab("Fraction") + coord_flip() +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values=annTypeColors, labels=rev(c("Up5kb", "Promoter", "5'UTR", "Exon", "Intron", "3'UTR", "Intergenic"))), nrow=1, rel_widths=c(0.3,0.9)),
  nrow=2, rel_heights=c(0.1,0.9)
)
gg_uniquelymissed


# Plot coverage around TSS
asdf  = genomation::annotateWithGeneParts(GR_results, geneRegionAnn)
asdf2 = genomation::getAssociationWithTSS(asdf)
GR_results2 = cbind(GR_results, asdf2[, c(2:3)])
GR_results2 = separate(GR_results2, col="Comparison", sep="\\.", into=c("assay", "algo"), remove=F)
GR_results3 = subset(GR_results2, dist.to.feature < 2500 & dist.to.feature > -2500)

gg_TSScov = ggplot(GR_results3, aes(x=dist.to.feature, y=depth, col=assay, fill=assay)) + 
  geom_hline(yintercept = 20, linetype="dashed", alpha=0.5) +
  stat_smooth(se=F, size=1, span=0.4, alpha=0.7) +
  theme_jon + theme(legend.position="right", legend.title=element_blank()) +
  ylab("Coverage") + xlab("Distance to TSS") +
  scale_color_manual(values = my.assay.cols) #+ facet_wrap(~assay, nrow=1)
gg_TSScov


# Plot coverage at CpG islands/shelves/shores
gg_CpGislandcov = ggplot(subset(resultsAnn, CpGAnnotation != "Intergenic"), aes(x=Comparison, y=depth)) +
  #geom_density(aes(fill=Comparison), alpha=0.7) +
  geom_violin(aes(fill=CpGAnnotation), trim = T, draw_quantiles = c(0.5), position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 20, linetype="dashed", alpha=0.7) +
  scale_fill_manual(values=c("CpG shelves"="#a1d99b", "CpG shores"="#74c476", "CpG islands"="#31a354", "Intergenic"="#bab3e3")) +
  #scale_color_manual(values = replicate(106, "#ebebeb")) +
  scale_x_discrete(labels = c("EMSeq"="EM", "MethylSeq"="MS", "SPLAT"="SP", "TruSeq"="TS", "TrueMethyl"="TM", "Nanopore"="NP")) +
  theme_jon + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5), legend.position="bottom", legend.title=element_blank(), legend.text = element_text(size=10)) +
  xlab("") + ylab("Mean Depth") + 
  coord_trans(ylim=c(0,50))
gg_CpGislandcov
