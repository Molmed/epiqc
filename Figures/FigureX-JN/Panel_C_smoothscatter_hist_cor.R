#Panels C-D of the Assay comparison figure:

library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(hexbin)
library(tidyr)
library(GGally)

# COLORS #
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

# LOAD DATA & ORGANIZE IT #
# FOR HG002, load the beta-values of all CpG sites >=10x coverage from down-sampled Bedgraph files to 20x. 
dat <- read.table(gzfile("/Users/jessicanordlund/Downloads/EPIQC/beta_HG002_LAB01_10x.txt.gz"), 
                 header= TRUE, fill=TRUE)  

indx <- grepl('EPIC', colnames(dat)) #Remove EPIC
dat <- dat[, !indx]

names(dat) <- c("TruSeq", "SPLAT", "EMSeq", "MethylSeq", "TrueMethyl") #Make clean names
count.cpgs.covered <- colSums(!is.na(dat)) #Count CPG sites with Beta-values


#Randomly sample 100k CGs: otherwise it takes too long
dat <- dat[sample(nrow(dat), 100000), ]

# format data for plotting
dat$CG <- row.names(dat) #add CG row name as ID 
tmp <- as_tibble(melt(dat, id="CG")) # Melt it 
tmp.cols <- my.assay.cols[match(names(dat), names(my.assay.cols))] # fix color order for plots
  
# Reduce DF- something happened with MethylSeq and TrueMethyl (all NA--- checking on this now).
df <- dat[,1:3]

#RMSE function 
RMSE = function(m, o){
  sqrt(mean((m - o)^2, na.rm=T))
}


# GGversion
p <- ggpairs(df, lower="blank", upper = "blank") 
seq <- 1:ncol(df)
  for (x in seq)
    for (y in seq) 
      if (y>x) 
        p <- putPlot(p, ggplot(df, aes_string(x=names(df)[x],y=names(df)[y])) + 
                       stat_density2d(aes(fill = ..density..^0.55), geom = "tile", contour = FALSE, n = 200) +
                       scale_fill_viridis_c(option = "plasma"), y,x)

for (x in seq) 
  for (y in seq) 
    if (x>y) {
      rmse.test <- round(RMSE(df[,x],df[,y]), 3)
      tmp.r <- round(cor(df[,x], df[,y], use = "pairwise.complete.obs", method="pearson"), 3)
      tmp.p <- round(cor(df[,x], df[,y], use = "pairwise.complete.obs", method="spearman"), 3)
      n.sites <- sum(complete.cases(df[,x],df[,y]))
      p <- putPlot(p, ggplot(df, aes_string(x=names(df)[x],y=names(df)[y])) + 
                     theme(panel.background = element_blank(), 
                           axis.title = element_blank(),
                           axis.text = element_blank(), 
                           axis.ticks = element_blank()) +
                     annotate("text", x=1, y=2,size =4, 
                              label= sprintf("r= %s \np= %s \nrmse= %s \n%s CpGs",tmp.r, tmp.p, rmse.test, n.sites)), y,x)}
p

