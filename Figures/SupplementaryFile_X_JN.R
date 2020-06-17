# Supplementary Figure for the Assay Comparison:
# Assay comparison by genome. 
# Based on Downsampled BedGraphs to 10x average coverage
# CpG site inclusion requirement is >= 5x coverage


library(data.table)
library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(hexbin)
library(tidyr)
library(GGally)

setwd("/Users/jessicanordlund/Downloads/EPIQC/for_manus")


# FUNCTIONS #
RMSE = function(m, o){
  sqrt(mean((m - o)^2, na.rm=T))
}

my.names <- data.frame(OLD= c("OXBS",
                              "EMSeq",
                              "MethylSeq",
                              "SPLAT",  
                              "TruSeq"), 
                       NEW =c("TrueMethyl", 
                              "EMSeq", 
                              "MethylSeq", 
                              "SPLAT", 
                              "TruSeq"))

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
# OBS TO DO, make a loop to parse through all genomes. 

# load the beta-values of all CpG sites >=10x coverage from down-sampled Bedgraph files to 20x. 
f <- '/Users/jessicanordlund/Downloads/EPIQC/10xDS/beta_HG001_LAB01_5x.txt.gz'
my.dat <- fread(f)

tmp.name <- sub("_LAB01_[[:digit:]]+x.txt.gz", "", basename(f))
tmp.name <- sub("beta_", "", tmp.name)

my.dat <- as.data.frame(my.dat)
indx <- grepl('EPIC', names(my.dat)) #Remove EPIC
my.dat<- my.dat[, !indx]

names(my.dat) <- gsub("_HG.*", "", names(my.dat))
names(my.dat) <- my.names$NEW[match(names(my.dat), my.names$OLD)] # clean names
my.dat <- my.dat[,order(names(my.dat))] #alphabetical order


# Calculate genome-wide correlations, RMSE and Pair-wise overlapping sites
count.cpgs.covered <- colSums(!is.na(my.dat)) #COunt CPG sites 
my.cor.p <- round(cor(my.dat, use = "pairwise.complete.obs", method="pearson"), 3)
my.cor.s <- round(cor(my.dat, use = "pairwise.complete.obs", method="spearman"), 3)

cor.p <- reshape2::melt(my.cor.p)
cor.s <- reshape2::melt(my.cor.s)

# Calc RMSE
combs <- combn(1:ncol(my.dat), 2) 
datalist = list()

for (i in 1:ncol(combs)) {
  x <- combs[1,i]
  y <- combs[2,i]
  datalist[[i]] <- c(names(my.dat)[x], names(my.dat)[y], 
                     round(RMSE(my.dat[,x],my.dat[,y]), 3),
                     sum(!is.na(my.dat[,x] & my.dat[,y])))
}

rmse_data = as.data.frame(do.call(rbind, datalist))
names(rmse_data) <- c("Lib1", "Lib2", "rmse", "nsites")


#Randomly sample 10k CGs to plot the smooth scatter: otherwise it takes too long to plot
df <- as.data.frame(my.dat[sample(nrow(my.dat), 10000), ])


# PLOT
p <- ggpairs(df, lower="blank", upper = "blank") 
seq <- 1:ncol(df)
for (x in seq)
  for (y in seq) 
    if (y>x) 
      p <- putPlot(p, ggplot(df, aes_string(x=names(df)[x],y=names(df)[y])) + 
                     stat_density2d(aes(fill = ..density..^0.55), geom = "tile", contour = FALSE, n = 200) +
                     theme(axis.text.x = element_text(angle = 90))+
                     scale_fill_viridis_c(option = "plasma"), y,x)

for (x in seq) 
  for (y in seq) 
    if (x>y) {
      rmse.test <- as.character(rmse_data$rmse[rmse_data$Lib1 %in% names(df)[x] & rmse_data$Lib2 %in% names(df)[y] |
                   rmse_data$Lib2 %in% names(df)[x] & rmse_data$Lib1 %in% names(df)[y]])
      n.sites <- as.character(rmse_data$nsites[rmse_data$Lib1 %in% names(df)[x] & rmse_data$Lib2 %in% names(df)[y] |
                  rmse_data$Lib2 %in% names(df)[x] & rmse_data$Lib1 %in% names(df)[y]])
      tmp.p <- cor.p$value[cor.p$Var1 %in% names(df)[x] & cor.p$Var2 %in% names(df)[y]]
      tmp.s <- cor.s$value[cor.s$Var1 %in% names(df)[x] & cor.s$Var2 %in% names(df)[y]]
     
        
      p <- putPlot(p, ggplot(df, aes_string(x=names(df)[x],y=names(df)[y])) + 
                     theme(panel.background = element_blank(), 
                           axis.title = element_blank(),
                           axis.text = element_blank(), 
                           axis.ticks = element_blank()) +
                     annotate("text", x=1, y=2,size =3, 
                              label= sprintf("r= %s \np= %s \nrmse= %s \n%s CpGs",tmp.p, tmp.s, rmse.test, n.sites)), y,x)}


png(filename= sprintf("%s_AssayComparison.png", tmp.name), 
    width=1300, height=850, res=200, bg = "transparent")
p
dev.off()


