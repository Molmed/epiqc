# Sequencing/microarray comparison plots
library(variancePartition)
library(data.table)
library(dplyr)
library(ggfan)
library(scales)
library(colorspace)
library(grid)
library(gridExtra)
library(gtable)
library(ggpubr)
library(viridis)
library(RColorBrewer)
theme_set(theme_bw(base_size = 9))
setwd('/path/to/epiqc/dir/')

###############################
##### Read in/format data #####
###############################
load('./results/VariancePartition_MicroarrayVsSequencing_Downsampled20x_BySeqPlatform.RData')
load('./results/VariancePartition_MicroarrayOnly_ByNormPipeline.RData')
load('./data/sequencing/BSseqData_Downsampled20x_MicroarraySitesOnly.RData')
load('./data/sequencing/MedianCoverage_SiteLevel_ByPlatform.RData')
load('./data/microarray/NormalizedDataSets.RData')
load('./data/microarray/rgset.RData')
samplesheet<-fread('./data/microarray/Samplesheet_AllLabs_EpiQCandControls.txt',data.table=F)
rownames(samplesheet)<-samplesheet$Basename

# Change 'PromethION' to 'Nanopore' to be consistent with other manuscript figures
PromethION.array.vp$SeqPlatform<-'Nanopore'

# Make large, merged dataset of sequencing/microarray comparisons
all.vp<-do.call(rbind,list(MethylSeq.array.vp,TruSeq.array.vp,EMSeq.array.vp,SPLAT.array.vp,TrueMethyl.array.vp,PromethION.array.vp))
common.cgs<-unique(all.vp$cgsite)
# Get median VEs for each platform/variable
all.vp %>% group_by(SeqPlatform) %>% summarise(Platform=median(VE.Platform),CellLine=median(VE.SampleName.GIAB),Residual=median(VE.Residuals))
# Nanopore noticably worse than other platforms

# Make large, merged dataset of microarray concordance results
objs<-objects(pattern='(.*)\\.(none|pbc|rcp|swan)\\.vp')
microarray.vp<-do.call(rbind,mget(objs))
# Keep only funnorm/RCP data individually (highest median)
objs<-objs[objs!='funnorm.rcp.vp']
rm(list=objs)
# Keep only sites present in all microarray analyses
cgtab<-table(microarray.vp$cgsite)
microarray.vp<-microarray.vp[microarray.vp$cgsite %in% names(cgtab[cgtab==max(cgtab)]),]
rm(cgtab)
microarray.vp$WithinArrayNorm<-factor(microarray.vp$WithinArrayNorm,levels=c('none','pbc','swan','rcp'),labels=c('None','PBC','SWAN','RCP'))
microarray.vp$BetweenArrayNorm<-factor(microarray.vp$BetweenArrayNorm,levels=c('none','pquantile','dasen','funnorm','enmix','sesame','gmqn'),labels=c('None','pQuantile','dasen','funnorm','ENmix','SeSAMe','GMQN'))


##### Combine all HG002 beta values #####
# From all 6 sequencing assays (1 merged replicate, downsampled bedgraph each)
bsseq.betas<-bsseq@assays@data@listData$Beta
bsseq.betas<-bsseq.betas[,grepl('HG002',colnames(bsseq.betas))]
# Remove EMSeq merged replicates, we want lab1/lab2 separated for the beta value plot
bsseq.betas<-bsseq.betas[,-which(colnames(bsseq.betas)=='EMSeq_HG002')]
bsseq.betas<-reshape2::melt(bsseq.betas)
colnames(bsseq.betas)<-c('cgsite','Sample','Beta')
bsseq.betas$cgsite<-as.character(bsseq.betas$cgsite)
bsseq.betas$Assay<-gsub('_HG002.*','',bsseq.betas$Sample)
bsseq.betas$Lab<-ifelse(grepl('LAB02',bsseq.betas$Sample),'Lab 2','Lab 1')
bsseq.betas$Replicate<-'Merged Replicates'

# All 3 microarray replicates (2 from Lab A, 1 from Lab B)
microarray.samples<-with(samplesheet,Basename[!is.na(SampleName.GIAB) & SampleName.GIAB=='HG002'])
microarray.betas<-reshape2::melt(funnorm.rcp.beta[,microarray.samples])
colnames(microarray.betas)<-c('cgsite','Sample','Beta')
microarray.betas$cgsite<-as.character(microarray.betas$cgsite)
microarray.betas$Assay<-'Microarray'
microarray.betas$Lab<-ifelse(microarray.betas$Sample %in% c('203175700016_R07C01','203175700087_R07C01'),'Lab A','Lab B')
microarray.betas$Replicate<-ifelse(microarray.betas$Sample %in% c('203175700016_R07C01','201959750090_R04C01'),'Replicate 1','Replicate 2')


# Combine into one dataframe w/ all HG002 information
hg002.all<-rbind(bsseq.betas,microarray.betas)
# Rename PromethION to Nanopore for consistency with other manuscript figures
hg002.all$Assay<-gsub('PromethION','Nanopore',hg002.all$Assay)
hg002.all$Label<-with(hg002.all,paste(Assay,Lab,Replicate,sep='\n'))
rm(bsseq.betas)
rm(microarray.betas)

##################################
##### Set plot colors/themes #####
##################################
# Color coding for all EpiQC plots by platform/sample
platform.colors<-c('#377eb8','#e41a1c','#4daf4a','#984ea3','#a65628','#4bcbde')
names(platform.colors)<-c('MethylSeq','TruSeq','EMSeq','SPLAT','TrueMethyl','Nanopore')
sample.colors<-c('#3aa142','#bdd7e7','#6baed6','#2171b5','#fcae91','#fb6a4a','#cb181d')
names(sample.colors)<-paste0('HG00',1:7)

# ggplot themes
theme_linedraw2 = theme_linedraw() +
  theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
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

##############################
##### Plotting functions #####
##############################
#### Individual CpG site plots (not included in manuscript) #####
MakeCGplot<-function(cgname) {
  seq.betavals<-data.frame(seq.beta=bsseq@assays@data@listData$Beta[cgname,],seq.coverage=bsseq@assays@data@listData$Coverage[cgname,])
  seq.betavals$SeqPlatform<-gsub('_.*','',rownames(seq.betavals))
  seq.betavals$SampleName.GIAB<-gsub('.*_','',rownames(seq.betavals))
  microarray.betavals<-data.frame(Basename=colnames(funnorm.rcp.beta),SampleName.GIAB=samplesheet[colnames(funnorm.rcp.beta),'SampleName.GIAB'],lab=samplesheet[colnames(funnorm.rcp.beta),'Lab'],microarray.beta=funnorm.rcp.beta[cgname,],stringsAsFactors=F)
  all.betavals<-merge(seq.betavals,microarray.betavals,by='SampleName.GIAB')
  all.betavals$SeqPlatform[all.betavals$SeqPlatform=='PromethION']<-'Nanopore'
  vp.res<-all.vp[all.vp$cgsite==cgname,]
  # Colored by cell line
  ggplot(data=all.betavals,aes(x=microarray.beta,y=seq.beta)) + stat_smooth(method='lm',alpha=0.2,colour='grey20') + geom_point(aes(size=seq.coverage,fill=SampleName.GIAB),alpha=0.8,shape=21) + scale_fill_manual(values=sample.colors,name='Cell line') + scale_size_continuous(name='Coverage') + facet_wrap(~SeqPlatform) + theme_jon + geom_label(data=vp.res,aes(x=-Inf,y=Inf,label=sprintf('VE by cell\nline: %.1f%%',VE.SampleName.GIAB*100)),hjust=-0.3,vjust=1.3)  + guides(fill = guide_legend(title.position="top", title.hjust = 0.5,nrow=2,override.aes=list(size=3)),size = guide_legend(title.position="top", title.hjust = 0.5,nrow=1)) + scale_x_continuous(name='Microarray beta value',labels=percent,limits=range(c(all.betavals$microarray.beta, all.betavals$seq.beta))) + scale_y_continuous(name='Sequencing beta value',labels=percent,limits=range(c(all.betavals$microarray.beta, all.betavals$seq.beta)))
}

#### Fan plots (supplementary figures comparing sequencing/microarray data) ####
MakeFanPlot<-function(platform=NULL,type='MicroarrayVar',probes='all') {
  if (type=='BetaValues') {
    if (probes=='high-varying') {
      plotdat<-hg002.plot[hg002.plot$cgsite %in% highvar.sites,]
      plotcolor<-'#66C2A5'
    } else if (probes=='low-varying') {
      plotdat<-hg002.plot[!(hg002.plot$cgsite %in% highvar.sites),]
      plotcolor<-'#FC8D62'
    } else if (probes=='all') {
      plotdat<-hg002.plot
      plotcolor<-'#66C2A5'
    }
  } else {
    plotdat<-get(sprintf('%s.array.vp',platform))
    plotdat<-plotdat[common.cgs,]
    if (platform=='PromethION') {
      plotcolor<-as.character(platform.colors['Nanopore'])
    } else {
      plotcolor<-as.character(platform.colors[platform])
    }
  }
  plot.intervals<-seq(0,0.9,by=0.1)
  if (type=='Coverage') {
    plotdat$CoverageCat<-round(plotdat$MedianCoverage)
    # After looking at the coverage histograms above for HG002,
    # we can see that we don't have enough observations to estimate
    # the plot intervals above coverage of 40
    lower.limit<-3
    upper.limit<-40
    summarized <- plotdat %>% data.frame() %>% filter(CoverageCat<upper.limit & CoverageCat>=lower.limit) %>% calc_quantiles(intervals=plot.intervals,x_var='CoverageCat',y_var='VE.SampleName.GIAB') %>% mutate(RowVar='Median coverage\n(sequencing data)',ColVar=ifelse(platform=='PromethION','Nanopore',platform))
  } else if (type=='MicroarrayVar') {
    plotdat$VarianceCat<-round(plotdat$MicroarrayVariance,3)
    # We have enough information to estimate plot quantiles for all sites
    # with microarray variance<0.116
    lower.limit<-0
    upper.limit<-0.116
    summarized <- plotdat %>% data.frame() %>% filter(VarianceCat<=upper.limit & VarianceCat>=lower.limit) %>% calc_quantiles(intervals=plot.intervals,x_var='VarianceCat',y_var='VE.SampleName.GIAB') %>% mutate(RowVar='CpG site variance\n(microarray data)',ColVar=ifelse(platform=='PromethION','Nanopore',platform))
  } else if (type=='SeqVar') {
    plotdat$VarianceCat<-round(plotdat$SeqVariance,3)
    # We have enough information to estimate plot quantiles for all sites
    # with variance<0.14 across all platforms (based on histogram above)
    lower.limit<-0
    upper.limit<-0.14
    summarized <- plotdat %>% data.frame() %>% filter(VarianceCat<=upper.limit & VarianceCat>=lower.limit) %>% calc_quantiles(intervals=plot.intervals,x_var='VarianceCat',y_var='VE.SampleName.GIAB') %>% mutate(RowVar='CpG site variance\n(sequencing data)',ColVar=ifelse(platform=='PromethION','Nanopore',platform))
  } else if (type=='BetaValues') {
    plotdat$Beta.x<-round(plotdat$Beta.x,2)
    x.freq<-table(plotdat$Beta.x)
    # Need a minimum of 20 observations to estimate quantiles
    plotdat<-plotdat[!(plotdat$Beta.x %in% as.numeric(names(x.freq[x.freq<20]))),]
    summarized <- plotdat %>% data.frame() %>% calc_quantiles(intervals=plot.intervals,x_var='Beta.x',y_var='Beta.y')
  } else {
    stop("type must be one of 'Coverage', 'MicroarrayVar', 'SeqVar', or 'BetaValues'")
  }
  # Plot data
  p<-ggplot(data=summarized,aes(x=x,y=y,quantile=quantile)) + geom_fan() + geom_interval(colour=plotcolor,intervals=c(0,0.5,0.9))
  # Add theme elements
  p <- p + scale_fill_gradient(low=lighten(plotcolor,amount=0.2),high='grey95',guide=F) + theme_jon + theme(legend.position='none', axis.text.x=element_text(angle=0),axis.title.x=element_blank(),axis.title.y=element_blank(),panel.grid.major=element_line(colour='grey70')) + scale_linetype_manual(values=c(1,2,3)) + scale_y_continuous(labels=percent,limits=c(0,1))
  if (type=='BetaValues') {
    p <- p + scale_x_continuous(labels=percent,limits=c(0,1))
  } else {
    p <- p + facet_grid(ColVar~RowVar)
  }
  return(p)
}


#######################
##### Make plots  #####
#######################
#### Individual CpG site plots (not included in manuscript) #####
# Example of poor TruSeq/array concordance, driven by low coverage
MakeCGplot('cg00003266')
# Example of poor Nanopore/array concordance, not driven by low coverage
MakeCGplot('cg00008893')
# Example of good concordance across the board; shows that results are largely driven
# by microarray data because they make up the majority of the samples (see Nanopore
# facet, which has high VE by cell line but clearly does not agree)
MakeCGplot('cg00001099')

# Plot Sequencing variance by platform to figure out where we have
# enough data to estimate the quantiles for the fan plot across
# all sequencing assays
ggplot(data=all.vp) + geom_hline(aes(yintercept=20),linetype=2,colour='grey50') + geom_histogram(aes(x=SeqVariance,colour=SeqPlatform, fill=SeqPlatform),binwidth=0.001,alpha=0.7) + facet_wrap(~SeqPlatform) + scale_fill_manual(values=platform.colors) + scale_colour_manual(values=platform.colors) + theme_jon + scale_y_log10() + theme(axis.text.x=element_text(angle=0)) + scale_x_continuous(breaks=10)

##### Figure 6 #####
### Figure 6a ###
sumstats <- microarray.vp %>% group_by(BetweenArrayNorm,WithinArrayNorm) %>% summarise(median.cellline=median(SampleName.HG),median.lab=median(Lab))  %>% ungroup()
# Manually add NAs for pipelines that don't work
sumstats <- sumstats %>% add_row(BetweenArrayNorm=c('SeSAMe','pQuantile'),WithinArrayNorm=c('SWAN','SWAN'),median.cellline=c(NA,NA),median.lab=c(NA,NA))

plotdat<-merge(microarray.vp,sumstats[,c('BetweenArrayNorm','WithinArrayNorm','median.cellline')],by=c('BetweenArrayNorm','WithinArrayNorm'))
# Order norm methods for plot (by publication date)
plotdat$BetweenArrayNorm<-factor(plotdat$BetweenArrayNorm,levels=c('None','pQuantile','dasen','funnorm','ENmix','SeSAMe','GMQN'))
plotdat$WithinArrayNorm<-factor(plotdat$WithinArrayNorm,levels=c('None','PBC','SWAN','RCP'))
sumstats$BetweenArrayNorm<-factor(sumstats$BetweenArrayNorm,levels=c('None','pQuantile','dasen','funnorm','ENmix','SeSAMe','GMQN'))
sumstats$WithinArrayNorm<-factor(sumstats$WithinArrayNorm,levels=c('None','PBC','SWAN','RCP'))
# No color or alpha scale
fig5a<-ggplot() + geom_density(data=plotdat, aes(x=SampleName.HG,fill=median.cellline)) + facet_grid(BetweenArrayNorm~WithinArrayNorm) + scale_x_continuous(label=percent,name='DNA methylation variance explained by cell line') + scale_fill_gradientn(colours=brewer.pal(9,'OrRd')[3:7],na.value=NA,name='Median') + geom_vline(data=sumstats,aes(xintercept=median.cellline),linetype=2) + geom_label(data=sumstats[!is.na(sumstats$median.cellline),],aes(label=sprintf('Median: %.0f%%',median.cellline*100),x=0.3,y=4),size=3,label.padding = unit(0.1, "lines")) + geom_text(data=sumstats[is.na(sumstats$median.cellline),],aes(x=0.5,y=3),label='Incompatible\nmethods',size=3) + theme_bw() + theme(panel.spacing.x = unit(0.5, "lines"), plot.title.position = "plot", strip.text=element_text(size=6))
ggsave(fig5a,filename='/projects/abv/GRC/lentsx/epiqc/figures/Figure5a.png',width=8,height=4,dpi=400)

### Figure 6b ###
# Determine microarray variance threshold
snp.betas<-getSnpBeta(rgset)
# Restrict to EpiQC samples (remove methylated/unmethylated controls)
epiqc.samples<-samplesheet$Basename[samplesheet$Type=='EpiQC']
snp.betas<-snp.betas[,epiqc.samples]
snpb.melt<-reshape2::melt(snp.betas)
snpb.melt$Var2<-as.character(snpb.melt$Var2)
snpb.melt<-merge(snpb.melt,samplesheet,by.x='Var2',by.y='Basename')

fig5b <- ggplot(data=snpb.melt) + geom_point(aes(y=Var1,x=value,colour=Lab),size=2,alpha=0.6) + coord_flip() + scale_colour_brewer(palette='Dark2',name='') + scale_y_discrete(name='') + scale_x_continuous(name='Raw beta value (no normalization)',label=percent) + theme_minimal() + theme(legend.position='bottom',axis.text.x=element_text(size=6,angle=90),axis.text.y=element_text(size=6), axis.title=element_text(size=12), legend.text=element_text(size=14),legend.margin=margin(-0.25,0,0,-1,'cm'), panel.grid.minor=element_blank()) + guides(colour = guide_legend(override.aes = list(size =4)))
ggsave(fig5b,filename='/projects/abv/GRC/lentsx/epiqc/figures/Figure5b.png',width=6,height=4,dpi=400)

# Call genotypes (naive method, fine because we have clean clusters)
snp.genotype<-matrix(NA,nrow=nrow(snp.betas),ncol=ncol(snp.betas),dimnames=list(rownames(snp.betas),colnames(snp.betas)))
snp.genotype[snp.betas<0.25]<-0
snp.genotype[snp.betas>0.75]<-2
snp.genotype[snp.betas>=0.25 & snp.betas<=0.75]<-1

rs.variances<-expand.grid(rownames(snp.betas),unique(samplesheet$Lab),c(0,1,2),stringsAsFactors=F)
colnames(rs.variances)<-c('SNP','Lab','Genotype')

# Calculate variance within each lab and genotype cluster
rs.variances$beta.var<-apply(rs.variances,1,function(x) {samples<-names(which(snp.genotype[x[[1]],colnames(snp.genotype) %in% samplesheet$Basename[samplesheet$Lab==x[[2]]]]==x[[3]])); geno.var<-ifelse(length(samples)==0,NA,var(snp.betas[x[[1]],samples]))})

summary(rs.variances$beta.var)
# Calculate 95th percentile to use as "low-varying" threshold
varthresh<-as.numeric(quantile(rs.variances$beta.var,na.rm=T,0.95))

# Calculate variances in funnorm/RCP normalized data for all probes
funnorm.rcp.vars<-apply(funnorm.rcp.beta[complete.cases(funnorm.rcp.beta),],1,var)

fig5c <- ggplot(data=rs.variances) + geom_histogram(aes(x=beta.var,fill=Lab),bins=50,alpha=0.8) +  scale_x_continuous(name='Methylation (raw beta value) variance') + scale_y_continuous(name='# of genotype clusters') + scale_fill_brewer(palette='Dark2',guide=F) + geom_vline(aes(xintercept=varthresh),linetype=2,colour='grey35',size=0.75) + annotate("segment", x=0.00155,xend=0.0013,y=45,yend=45,arrow=arrow(length = unit(0.1, "cm")), color="grey35") + annotate("text", x=0.0016,y=45, label="95th percentile", color = "grey35",size=2.75, hjust=0) + theme_bw() + theme(axis.title=element_text(size=8),axis.text=element_text(size=6))
ggsave(fig5c,filename='/projects/abv/GRC/lentsx/epiqc/figures/Figure5c.png',width=4,height=2,dpi=400)


fig5d <- ggplot() + geom_histogram(aes(x=funnorm.rcp.vars),bins=60,fill='#E6AB02',colour='grey35',alpha=0.7,size=0.25) + geom_vline(aes(xintercept=varthresh),size=0.5,linetype=5,colour='grey35') + annotate("rect", xmin = -Inf, xmax = max(rs.variances$beta.var,na.rm=T), ymin = 0, ymax = Inf, fill = "grey35", alpha = 0.4)  + annotate("segment", x=0.028,xend=0.005,y=325000,yend=325000,arrow=arrow(length = unit(0.1, "cm")), color="grey35") + annotate("text", x=0.03,y=325000, label="SNP probe 95th percentile", color = "grey35",size=3, hjust=0) + geom_segment(aes(x=0.04,xend=0.11,y=170000,yend=170000),arrow = arrow(length = unit(0.1, "cm"))) + geom_text(aes(x=0.02,y=200000,label='High-varying CpG sites'),size=3,hjust=0) + scale_x_continuous(name='Methylation variance (funnorm/RCP beta values)',expand=c(0.04,0))  + scale_y_continuous(name='# of CpG sites',labels=comma) + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=6))
ggsave(fig5d,filename='/projects/abv/GRC/lentsx/epiqc/figures/Figure5d.png',width=3.75,height=2,dpi=400)

funnorm.rcp.vp$MicroarrayVariance<-funnorm.rcp.vars[funnorm.rcp.vp$cgsite]
funnorm.rcp.vp$MicroarrayVarCat<-ifelse(funnorm.rcp.vp$MicroarrayVariance>=varthresh,'High-varying sites','Low-varying sites')
funnorm.rcp.vp$MicroarrayVarCat<-factor(funnorm.rcp.vp$MicroarrayVarCat,levels=c('Low-varying sites','High-varying sites'),labels=c(sprintf('Low-varying sites (N=%s)',comma(sum(funnorm.rcp.vp$MicroarrayVarCat=='Low-varying sites'))),sprintf('High-varying sites (N=%s)',comma(sum(funnorm.rcp.vp$MicroarrayVarCat=='High-varying sites')))))
#ggplot(data=funnorm.rcp.vp) + geom_histogram(aes(x=SampleName.HG),colour='grey35',fill='#ffff33',bins=40) + facet_wrap(~MicroarrayVarCat) + scale_x_continuous(name='% variance explained by cell line') + scale_y_continuous(name='# of CpG sites',labels=comma)
fig5e <- ggplot(data=funnorm.rcp.vp) + geom_histogram(aes(x=SampleName.HG),bins=50,fill='#E6AB02',colour='grey35',size=0.25) + facet_wrap(~MicroarrayVarCat) + scale_x_continuous(name='Variance explained by cell line (funnorm/RCP normalized data)',labels=percent) + scale_y_continuous(name='# of CpG sites',labels=comma) + theme_bw() + theme(strip.text=element_text(size=8),axis.title=element_text(size=9))
ggsave(fig5e,filename='/projects/abv/GRC/lentsx/epiqc/figures/Figure5e.png',width=6,height=2,dpi=400)


### Make fan plot of variance partition results (Figure 6b) ###
# Making the pre-calculated quantiles work with facet_grid was
# a pain, so we're just looping through and determining whether
# the facet label should be on the plot based on the position
plot.rows<-c('EMSeq','MethylSeq','PromethION','SPLAT','TrueMethyl','TruSeq')
plot.columns<-c('MicroarrayVar','SeqVar','Coverage')
fanplot.vp<-list()
for (r in plot.rows) {
  for (c in plot.columns) {
    p<-MakeFanPlot(platform=r,type=c)
    if (r!=plot.rows[[length(plot.rows)]]) {
      p <- p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
    }
    if (r!=plot.rows[[1]]) {
      p <- p + theme(strip.text.x = element_blank(), strip.background.x=element_blank())
    }
    if (c!=plot.columns[[1]]) {
      p <- p + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
    }
    if (c!=plot.columns[[length(plot.columns)]]) {
      p <- p + theme(strip.text.y = element_blank(), strip.background.y=element_blank())
    }
    fanplot.vp[[length(fanplot.vp)+1]]<-p
  }
}
# Make shared legend
legplot<-ggplot(data=data.frame(x=rep(1:10,times=3),y=rep(c(1,2,3),each=10),line=rep(c('Median','0.5','0.9'),each=10))) + geom_line(aes(x=x,y=y,linetype=line)) + scale_linetype_manual(values=c(1,2,3),breaks=c('Median','0.5','0.9'),name='Interval') + theme_jon + theme(legend.position='bottom')
legend = gtable_filter(ggplotGrob(legplot), "guide-box")

# Make full figure
ggsave(arrangeGrob(arrangeGrob(grobs=fanplot.vp,ncol=3,left=textGrob('Variance explained by cell line in merged sequencing/microarray dataset',rot=90,vjust=1),widths=c(2.2,1.9,2.1),heights=c(2.9,2,2,2,2,2.3)),legend,ncol=1,heights=c(20,1)),width=7,height=7,filename='/projects/abv/GRC/lentsx/epiqc/figures/CrossPlatformConcordance_ByVarianceCoverage_FanPlot.png',dpi=400)

#### Beta value comparison, by microarray variance (high or low) ####
# Again, we're plotting based on the position here; the lower left
# will be low-varying sites and upper right high-varying; the diagonal
# indicates the assay/replicate

plot.samples<-unique(hg002.all$Label)[1:3]
fanplot.beta<-list()
for (r in 1:(length(plot.samples)+1)) {
  for (c in 1:length(plot.samples)) {
    if(r==c) {
      p<-ggplot(data=hg002.all[hg002.all$Label==plot.samples[c],]) + geom_histogram(aes(x=Beta),fill='grey80',colour='grey30') + scale_y_continouos() + theme_jon + theme(axis.title=element_blank(),panel.grid.major=element_line(colour='grey70'))
    } else if (r==(c+1)) {
        p<-textGrob(plot.samples[r],gp=gpar(fontsize=7))
    } else if (r!=(c+1)) {
      hg002.plot<-merge(hg002.all[hg002.all$Sample==plot.samples[c],],hg002.all[hg002.all$Sample==plot.samples[r],],by=c('cgsite'))
      hg002.plot<-hg002.plot[hg002.plot$cgsite %in% common.cgs,]
      if(r<c) {
        #p<-MakeFanPlot(type='BetaValues',probes='high-varying')
        p<-ggplot() + theme_void()
      } else if (r>c) {
        #p<-MakeFanPlot(type='BetaValues',highvarying='low-varying')
        p<-MakeFanPlot(type='BetaValues',probes='all') + theme(plot.margin = unit(c(2,2,2,2),'pt'))
      }
      if (r!=length(plot.samples) & r!=1) {
        p<-p + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
      } else if (r==1) {
        p <- p + scale_x_continuous(position='top',labels=percent) + theme(axis.text.x=element_text(size=5))
      } else if (r==length(plot.samples)) {
        p <- p + scale_x_continuous(labels=percent) + theme(axis.text.x=element_text(size=5))
      }
      if (c!=length(plot.samples) & c!=1) {
        p <- p + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
      } else if (c==length(plot.samples)) {
        p <- p + scale_y_continuous(position='right',labels=percent) + theme(axis.text.y=element_text(size=5))
      } else if (c==1) {
        p <- p + scale_y_continuous(labels=percent) + theme(axis.text.y=element_text(size=5))
      }
    }
    fanplot.beta[[length(fanplot.beta)+1]]<-p
  }
}
# grid.arrange(arrangeGrob(grobs=fanplot.vp,ncol=3,left=textGrob('Variance explained by cell line in merged sequencing/microarray dataset',rot=90,vjust=1),heights=c(2.7,2,2,2,2,2)),legend,ncol=1,heights=c(20,1))
h.vec<-c(rep(1,length(plot.samples)),1.15)
w.vec<-c(1.25,rep(1,length(plot.samples)-1))
ggsave(arrangeGrob(arrangeGrob(grobs=fanplot.beta,heights=h.vec,widths=w.vec),legend,ncol=1,heights=c(25,1)),width=9,height=9,filename='/projects/abv/GRC/lentsx/epiqc/figures/HG002_BetaValue_Comparison.png',dpi=400)


#### Ternary plots ####
# There are conflicts between ggplot themes and ggtern themes, easiest to
# just load ggtern at the end and make this plot last
# This is part of figure 6, but I plotted it at the end because of the
# above conflicts
library(ggtern)
plotdat<-all.vp[,c('VE.SampleName.GIAB','VE.Platform','VE.Residuals','SeqPlatform')]
plotdat[,1:3]<-100*plotdat[,1:3]
fig6a<-ggtern(data=plotdat,aes(x=VE.SampleName.GIAB,y=VE.Platform,z=VE.Residuals),aes(x,y,z)) +
  stat_density_tern(geom="polygon",
                    n=100,#h=0.5,
                    base='identity',
                    aes(alpha = ..level..,fill=..level..),
                    bins=50,
                    na.rm = TRUE) +
  scale_fill_viridis(alpha=0.8,name='Level') +
  labs(x='',xarrow='VE by cell line (%)',y='',yarrow='VE by assay (%)',z='',zarrow='\nResidual variation (%)') +
  facet_wrap(~SeqPlatform,ncol=3) +
  theme_showarrows() + theme(panel.background = element_rect(fill = NA), panel.grid.major=element_line(size=0.2,linetype='dashed',colour='grey80'), panel.ontop=TRUE, panel.grid.minor=element_line(colour=NA),axis.text=element_text(vjust=0),legend.position='right') + scale_alpha_continuous() + guides(alpha=F)
ggsave(fig5f,filename='/projects/abv/GRC/lentsx/epiqc/figures/Figure6a.png',width=10,height=8,dpi=400)
