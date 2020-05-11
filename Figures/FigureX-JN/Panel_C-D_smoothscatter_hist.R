#Panels C-D of the Assay comparison figure:

library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)

my.assay.cols <- c("MethylSeq" = "#377eb8", 
                   "TruSeq" = "#e41a1c",
                   "EMSeq" = "#4daf4a",
                   "SPLAT" = "#984ea3",
                   "EPIC" = "#ff7f00",
                   "Arrays" = "#ffff33",
                   "TrueMethyl" = "#a65628")

my.genome.cols <- c("HG001" = "#3aa142",
                    "HG002" = "#bdd7e7",
                    "HG003" = "#6baed6",
                    "HG004" = "#2171b5",
                    "HG005" = "#fcae91",
                    "HG006" = "#fb6a4a",
                    "HG007" = "#cb181d")

# FOR HG002, load the beta-values of all CpG sites >=10x coverage from down-sampled Bedgraph files to 20x. 
dat <- read.table(gzfile("/Users/jessicanordlund/Downloads/EPIQC/beta_HG002_LAB01_10x.txt.gz"), 
                 header= TRUE, fill=TRUE)  
count.cpgs.covered <- colSums(is.na(dat))

indx <- grepl('EPIC', colnames(dat)) #Remove EPIC
dat <- dat[, !indx]

names(dat) <- c("TruSeq", "SPLAT", "EMSeq", "MethylSeq", "TrueMethyl") #Make clean names
count.cpgs.covered <- colSums(!is.na(dat))


#Randomly sample 1M CGs: otherwise it takes too long
dat <- dat[sample(nrow(dat), 1000000), ]

# format data for plotting
dat$CG <- row.names(dat) #add CG row name as ID 
tmp <- as_tibble(melt(dat, id="CG")) # Melt it 
tmp.cols <- my.assay.cols[match(names(dat), names(my.assay.cols))] # fix color order for plots
  
# Get data, calculate correlations
df <- dat[,1:3]
my.cor.p <- round(cor(df, use = "pairwise.complete.obs", method="pearson"), 2)
my.cor.s <- round(cor(df, use = "pairwise.complete.obs", method="spearman"), 2)

# ----------------------------------------------------#
# PLOT THE SMOOTH SCATTER PLOTS- PANEL C
# Build colors
buylrd = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
           "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026") 
myColRamp = colorRampPalette(c(buylrd))

# Get all info organized
combs <- combn(1:ncol(df), 2) 
my.cor.p<- melt(my.cor.p)
my.cor.s <- melt(my.cor.s)

#Make the plot
par(mfrow=c(2,2), mar=c(3,4,1,1), mgp=c(1.5, .5, 0))

for(i in 1:ncol(combs)){
  x.name <-names(df)[combs[1,i]]
  y.name <- names(df)[combs[2,i]]
  cor.text.nr <- which(my.cor.p$Var1 %in% x.name & my.cor.p$Var2 %in% y.name)
  legend.text <- c(sprintf("r=%s",my.cor.p$value[cor.text.nr]), sprintf("p=%s",my.cor.s$value[cor.text.nr]))
  
  smoothScatter(x=df[,combs[1,i]], y=df[,combs[2,i]],
                xlab=x.name,ylab=y.name,
                colramp=myColRamp)
  
  legend("topleft", legend=legend.text, cex=0.9, bty="n", text.col="white")
  }

#-------------------------------------------------------------------------------#



#
# PLOT the Beta-value distribution: Panel D
tmp <- rename(tmp, Assay = variable)

ggplot(tmp, aes(x=value, color=Assay))+
  geom_density() +
  scale_color_manual(values= tmp.cols, 
                    labels = paste0(levels(tmp$Assay), " (", round(count.cpgs.covered/1000000, 1), ")"))+
  xlab("beta-value") +
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


  

