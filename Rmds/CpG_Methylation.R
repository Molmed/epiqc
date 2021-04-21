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
library(GenomicRanges)
theme_linedraw2 = theme_linedraw() +  theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
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

splitGet=function (strings, sep, n)
{
  return ( sapply(strsplit(strings, sep), "[[", n) )
}

read.tsv = function(file, header=T, ...)
{
  return(read.table(file, sep="\t", header=header, stringsAsFactors=F,  na.strings = c("NA","#N/A","N/A"), ...))
}


read.tsv2 = function(file, header=T, row.names=F, ...)
{
  if(row.names==F)
  {
    library(data.table)
    df = fread(file, stringsAsFactors=F, na.strings = c("NA","#N/A","N/A"), data.table=F, ...)
    colnames(df) = make.names(colnames(df))
    return(df)
  }else
  {
    return(read.table(file, sep="\t", header=header, stringsAsFactors=F,  na.strings = c("NA","#N/A","N/A"), ...))
  }
}

setRemoveRownames = function(df, rownameColumn=1)
{
  rownames(df) = df[,rownameColumn]
  df = df[,-rownameColumn, drop=F]
  return(df)
}


getOverlapCpGDF = function(peakGR, cpgGR, shift = 0, minHitsPerPeak = 0)
{
  ov = findOverlaps(peakGR, cpgGR)
  dfM = data.frame( peakID = peakGR$name[queryHits(ov)],
                    peakGene = peakGR$gene[queryHits(ov)],
                    chr = peakGR@seqnames[queryHits(ov)],
                    peakStart = peakGR@ranges@start[queryHits(ov)],
                    peakEnd = peakGR@ranges@start[queryHits(ov)]+peakGR@ranges@width[queryHits(ov)]-1,
                    peakWidth = peakGR@ranges@width[queryHits(ov)],
                    cpgId = cpgGR$name[subjectHits(ov)],
                    cpgStart = cpgGR@ranges@start[subjectHits(ov)])
  if(minHitsPerPeak > 0)
  {
    peakCounts = table(dfM$peakID)
    dfM = dfM[dfM$peakID %in% names(peakCounts)[peakCounts >= minHitsPerPeak], ]
  }
  
  
  dfM$peakMid = dfM$peakStart + (dfM$peakWidth/2)
  dfM$cpgRelative = (dfM$cpgStart - dfM$peakMid) / (dfM$peakWidth/2)
  dfM$cpgRelative = dfM$cpgRelative + shift
  return(dfM)
}

getGeneRegionCpG = function(cpgGR, peaksList, minHitsPerPeak = 3, binScaling = 2)
{
  require(GenomicRanges)
  require(reshape2)
  require(dplyr)
  
  dfAll = data.frame()
  binScalingVals = if( length(binScaling) != length(peaksList)) rep(binScaling[1], length(peaksList)) else binScaling
  
  for(i in 1:length(peaksList))
  {
    peaks = peaksList[[i]]
    peaks$name = paste(peaks$chr, peaks$start, peaks$end, sep=".")
    peakType = names(peaksList)[i]
    peakGR = GRanges(peaks$chr, IRanges(start=peaks$start, end=peaks$end), name=peaks$name, gene=peaks$gene)
    
    dfOverlap = getOverlapCpGDF(peakGR, cpgGR, minHitsPerPeak = minHitsPerPeak, shift=(1+(i-1)*2))
    if(nrow(dfOverlap) == 0)
      next
    dfOverlap$Type = peakType
    dfOverlap$cpgRelative2 = round(dfOverlap$cpgRelative*binScalingVals[i], 1) / binScalingVals[i]
    dfAll = rbind(dfAll, dfOverlap)
  }
  return(dfAll)
}


rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
r  = rf(32)

# set up dataframe
df = read.csv("cormatrix/corMatrix_METH_ds20.chr1.csv.gz")
df = separate(df, col="locus", sep="\\.", into=c("chr","start","end"), remove=F)
df$chrStart = paste0(df$chr, ".", df$start)
df = df[, -c(2:4)]
dfm = melt(df)
colnames(dfm) = c("locus", "chrStart", "assay", "meth")
dfm$assay = factor(dfm$assay, levels=c("EMSeq", "MethylSeq", "Nanopore", "TrueMethyl", "SPLAT", "TruSeq"))

# Methylation Distribution Plot
CpG_meth = gghistogram(dfm, x="meth", color="assay", fill="assay", 
            palette=my.assay.cols, alpha=0.7, 
            xlim=c(0,100), bins=20, binwidth = 5, 
            xlab="Methylation %", ylab="Frequency",
            legend=0,
            facet.by = "assay", nrow=2,
            add="median",
            ggtheme=theme_bw()) +
  theme_jon + theme(legend.position="none", strip.text.y = element_text(angle=0, size=10))
CpG_meth

# Correlation plots
correlation = cor(df[2:7], use="complete.obs")
correlation_m = melt(correlation)
correlation_m$value = round(correlation_m$value, 2)
gg_20x_corr = ggplot(correlation_m, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  geom_text(aes(label = value), size=4, color="#000000") +
  theme_jon + theme(legend.position="none") +
  #scale_fill_viridis(option = "plasma", ) +
  scale_fill_viridis(option = "plasma", limits = c(0.5, 1), oob = scales::squish) +
  labs(fill = "Pearson") +
  xlab("") + ylab("") + ggtitle("All CpGs")
gg_20x_corr


df$mean = rowMeans(df[, 2:7])
correlation20to80 = cor(subset(df, mean < 80 & mean > 20)[2:7], use="complete.obs")
correlation20to80_m = melt(correlation20to80)
correlation20to80_m$value = round(correlation20to80_m$value, 2)
gg_20x_corr20to80 = ggplot(correlation20to80_m, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  geom_text(aes(label = value), size=4, color="#000000") +
  theme_jon + theme(legend.position="right") +
  scale_fill_viridis(option = "plasma", limits = c(0.5, 1), oob = scales::squish) +
  labs(fill = "Pearson") +
  xlab("") + ylab("") + ggtitle("20-80% Methylated CpGs")
gg_20x_corr20to80


# mBias plot
df_mbias = read.csv("tables/mBias_all.csv")
df_mbias$pct = df_mbias$nMethylated / (df_mbias$nMethylated + df_mbias$nUnmethylated) * 100
df_mbias = subset(df_mbias, strand %in% c("OT", "OB"))
df_mbias$strand = factor(df_mbias$strand, levels=c("OT", "OB"))
df_mbias = separate(df_mbias, col="sample", into=c("assay","genome","lab","rep"), remove=F)
df_mbias$assayStrand = paste0(df_mbias$assay, "_", df_mbias$strand)
df_mbias$assay = factor(df_mbias$assay, levels=c("EMSeq", "MethylSeq", "TrueMethylBS", "TrueMethylOX", "SPLAT", "TruSeq"))

mbias = ggplot(subset(df_mbias, algo=="bwameth"), aes(position, pct)) +
  geom_line(stat="smooth", aes(color=assay), se = F, span = 0.7, alpha=0.8, size=1.5) +
  geom_vline(xintercept = c(20,130), linetype="dashed", alpha=0.7) +
  scale_x_continuous(breaks=c(0, 20, 50, 100, 130, 150)) +
  scale_color_manual(values = my.assay.cols) +
  theme_jon + theme(legend.title = element_blank()) +
  coord_trans(ylim=c(55,75)) +
  facet_grid(strand~read, labeller = labeller(strand = c("OT"="Orig Top", "OB"="Orig Bottom"), read = c("1"= "Read 1", "2"= "Read 2"))) +
  xlab("Base") + ylab("% Methylation")
mbias

# ONT vs WMS
correlation = cor(df$Nanopore, df$EMSeq, use="complete.obs")
plotTitle = sprintf(paste0("\nr = %.3f"), correlation)
gg_PvsEM = ggMarginal(ggplot(df, aes(y=Nanopore, x=EMSeq)) + geom_bin2d(bins=25, show.legend=F) + geom_point(col="transparent") + scale_fill_gradientn(colors=r, trans="log10") + xlab("EMSeq") + ylab("Nanopore") + theme_jon + theme(plot.title=element_text(size=13), axis.text.y=element_text(size=13), axis.text.x=element_text(size=13), axis.title.x=element_text(size=13), axis.title.y=element_text(size=13)) + ggtitle(plotTitle), type="histogram", margins = "x") 
correlation = cor(df$Nanopore, df$MethylSeq, use="complete.obs")
plotTitle = sprintf(paste0("\nr = %.3f"), correlation)
gg_PvsMS = ggMarginal(ggplot(df, aes(y=Nanopore, x=MethylSeq)) + geom_bin2d(bins=25, show.legend=F) + geom_point(col="transparent") + scale_fill_gradientn(colors=r, trans="log10") + xlab("MethylSeq") + ylab("") + theme_jon + theme(plot.title=element_text(size=13), axis.text.y=element_text(size=13), axis.text.x=element_text(size=13), axis.title.x=element_text(size=13), axis.title.y=element_text(size=13)) + ggtitle(plotTitle), type="histogram", margins = "x") 
correlation = cor(df$Nanopore, df$SPLAT, use="complete.obs")
plotTitle = sprintf(paste0("\nr = %.3f"), correlation)
gg_PvsSP = ggMarginal(ggplot(df, aes(y=Nanopore, x=SPLAT)) + geom_bin2d(bins=25, show.legend=F) + geom_point(col="transparent") + scale_fill_gradientn(colors=r, trans="log10") + xlab("SPLAT") + ylab("") + theme_jon + theme(plot.title=element_text(size=13), axis.text.y=element_text(size=13), axis.text.x=element_text(size=13), axis.title.x=element_text(size=13), axis.title.y=element_text(size=13)) + ggtitle(plotTitle), type="histogram", margins = "x") 
correlation = cor(df$Nanopore, df$TrueMethyl, use="complete.obs")
plotTitle = sprintf(paste0("\nr = %.3f"), correlation)
gg_PvsTM = ggMarginal(ggplot(df, aes(y=Nanopore, x=TrueMethyl)) + geom_bin2d(bins=25, show.legend=F) + geom_point(col="transparent") + scale_fill_gradientn(colors=r, trans="log10") + xlab("TrueMethyl") + ylab("") + theme_jon + theme(plot.title=element_text(size=13), axis.text.y=element_text(size=13), axis.text.x=element_text(size=13), axis.title.x=element_text(size=13), axis.title.y=element_text(size=13)) + ggtitle(plotTitle), type="histogram", margins = "x") 
correlation = cor(df$Nanopore, df$TruSeq, use="complete.obs")
plotTitle = sprintf(paste0("\nr = %.3f"), correlation)
gg_PvsTS = ggMarginal(ggplot(df, aes(y=Nanopore, x=TruSeq)) + geom_bin2d(bins=25, show.legend=F) + geom_point(col="transparent") + scale_fill_gradientn(colors=r, trans="log10") + xlab("TruSeq") + ylab("") + theme_jon + theme(plot.title=element_text(size=13), axis.text.y=element_text(size=13), axis.text.x=element_text(size=13), axis.title.x=element_text(size=13), axis.title.y=element_text(size=13)) + ggtitle(plotTitle), type="histogram", margins = "both")
ggPvsAll = plot_grid(gg_PvsEM, gg_PvsMS, gg_PvsSP, gg_PvsTM, gg_PvsTS, nrow=1)
ggPvsAll



# Metagene methylation plot

### Create GenomicRanges object
cpg = data.frame(chr=splitGet(df$locus, "\\.", 1), start=splitGet(df$locus, "\\.", 2))
cpg$start = as.numeric(as.character(cpg$start))
cpgGR = with(cpg, GRanges(chr, IRanges(start=start, width=1), name=paste0(chr, ".", start)))
### Get regional table (for hg38)
region = read.tsv2("hg38_withHeader_ucsc2016_allTableRegulatory.txt", sep='\t', header=TRUE)
region[1:5, 1:5]
### Create Gene Regions (and scaling for plot)
geneRegions = list()
geneRegions[["Up1k"]] = with(region[region$Type2=="up1k", ], data.frame(chr=chr, start=start, end=end, gene=geneSymbol))
geneRegions[["Promoter"]] = with(region[region$Type2=="regulatory", ], data.frame(chr=chr, start=start, end=end, gene=geneSymbol))
geneRegions[["5UTR"]] = with(region[region$Type2=="utr5", ], data.frame(chr=chr, start=start, end=end, gene=geneSymbol))
geneRegions[["Exon"]] = with(region[region$Type2=="exon", ], data.frame(chr=chr, start=start, end=end, gene=geneSymbol))
geneRegions[["Intron"]] = with(region[region$Type2=="intron", ], data.frame(chr=chr, start=start, end=end, gene=geneSymbol))
geneRegions[["3UTR"]] = with(region[region$Type2=="utr3", ], data.frame(chr=chr, start=start, end=end, gene=geneSymbol))
geneRegions[["Down5k"]] = with(region[region$Type2=="down5k", ], data.frame(chr=chr, start=start, end=end, gene=geneSymbol))
genePlotBreaks = seq(0,12,by=2)
genePlotLabels = c("Up 1kb", "Promoter", "5'UTR", "Exons", "Introns", "3'UTR", "Down 5kb")
genePlotBinScaling = c(3, 2, 2, 2, 2, 2, 3)
### Combine Gene Regions with CpG GR
cpgPos = getGeneRegionCpG(cpgGR, geneRegions, minHitsPerPeak = 3, binScaling=genePlotBinScaling)
cpgPos
### uniquify the relative position of CpGs and their associated genes
absMean = function(x) return(mean(abs(x), na.rm=T))
mean2 = function(x) return(mean(x, na.rm=T))
cpgPos2 = unique(cpgPos[, c("cpgId", "cpgRelative2", "peakGene")])
### merge with original CpG matrix
allCpG2 = merge(cpgPos2, df, by.x="cpgId", by.y="chrStart")
### summarise value at relative position along plot
allCpG3 = allCpG2[, -c(1,4)] %>%
  group_by(cpgRelative2) %>%
  summarise_each(funs(mean2), -peakGene)
### combine all samples into one column rather than a column for each sample
melted_allCpG3 = melt(allCpG3, "cpgRelative2")
### Plot!
ggplot(subset(melted_allCpG3, cpgRelative2>2), aes(x=cpgRelative2, y=value, col=variable, fill=variable)) + 
  geom_line(size=0.5) +
  geom_vline(xintercept=c(2,4,6,8,10,12,14), linetype = "longdash", alpha=0.7) + 
  theme_jon + theme(axis.text.x=element_text(size=10), legend.position="bottom", legend.title=element_blank()) +
  scale_color_manual(values=my.assay.cols) +
  ylab("Mean Methylation") + xlab("") +
  scale_x_continuous(breaks=c(3,4,5,7,9,11,13), labels=c("Promoter\n(Up 1kb)\n", "TSS", "5'UTR", "Exons", "Introns", "3'UTR", "Down\n5kb"))

### Plot TSS Methylation %
### this requires GR_results3 object from CpG_coverage.R
ggplot(GR_results3, aes(x=dist.to.feature, y=meth, col=assay, fill=assay)) + 
  stat_smooth(se=F, size=1, span=0.4, alpha=0.7, show.legend = F) +
  theme_jon + theme(legend.position="bottom", legend.title=element_blank()) +
  ylab("Methylation %") + xlab("Distance to TSS") +
  scale_color_manual(values = my.assay.cols)
