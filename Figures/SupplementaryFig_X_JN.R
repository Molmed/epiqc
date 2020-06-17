# Supplementary Figure for the Assay Comparison:
# Assay comparison by genome. 
# Based on Downsampled BedGraphs to 10x average coverage
# CpG site inclusion requirement is >= 5x coverage

suppressPackageStartupMessages()
library(data.table)
library(ggplot2) 
library(dplyr) 
library(reshape2)
library(RColorBrewer)
#library(hexbin)
library(tidyr)
library(GGally)
library(ggbeeswarm)
library(corrplot)
options(stringsAsFactors = FALSE)


# Set working directory
setwd("/proj/uppstore2017167/EPIQC/plots_for_manuscript/plot")


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

#----------------------------#
# GET THE DATA and format it #
#----------------------------#
dat.dir <- '/proj/uppstore2017167/EPIQC/plots_for_manuscript/10xDS_bedgraphs'
my.files <- list.files(dat.dir, recursive=F, full.names= T, pattern= "[[:digit:]]+x.txt.gz")

my.stats <- as.list(NULL)
my.cpgs <- as.list(NULL)
my.cor.p <- as.list(NULL)
my.cor.s <- as.list(NULL)
my.rmse <- as.list(NULL)
num.sites <- "100000"
large.df <- data.frame()[1:num.sites, ]



#----------------------------#
# Loop though by genome #
#----------------------------#

for (f in my.files) {
  my.dat <- fread(f)
  
  # 1 fix names & order
  {
  tmp.name <- sub("_LAB01_[[:digit:]]+x.txt.gz", "", basename(f))
  tmp.name <- sub("beta_", "", tmp.name)
  
  print(c(tmp.name, "is loaded"))
  
  my.dat <- as.data.frame(my.dat)
  indx <- grepl('EPIC', names(my.dat)) #Remove EPIC
  my.dat<- my.dat[, !indx]
  
  names(my.dat) <- gsub("_HG.*", "", names(my.dat))
  names(my.dat) <- my.names$NEW[match(names(my.dat), my.names$OLD)] # clean names
  my.dat <- my.dat[,order(names(my.dat))] #alphabetical order
  }

  # 2 add a subset to the large dataframe - for the corr plot
  {
  tmp <- my.dat[1:num.sites ,]
  names(tmp) <- paste(tmp.name,names(tmp), sep="_")
  large.df <- cbind(large.df, tmp)
  rm(tmp)
  }

  
  # 3 basic statistics
  {
  print("getting basic stats")
  my.stats[[tmp.name]] <- apply(my.dat, 2, summary, digits = 3) 
  my.cpgs[[tmp.name]] <- colSums(!is.na(my.dat))
  }
  
  
  # 4 calculating correlations
  {
  print("Calculating Correlation")
  cor.p <- round(cor(my.dat, use = "pairwise.complete.obs", method="pearson"), 3)
  cor.s <- round(cor(my.dat, use = "pairwise.complete.obs", method="spearman"), 3)
  
  my.cor.p[[tmp.name]] <- reshape2::melt(cor.p)
  my.cor.s[[tmp.name]] <- reshape2::melt(cor.s)
  }
  
  
  # 5 calculating RMSE 
  {
    print("Calculating RMSE")
    combs <- combn(1:ncol(my.dat), 2) 
    datalist = list()
  
    for (i in 1:ncol(combs)) {
      x <- names(my.dat)[combs[1,i]]
      y <- names(my.dat)[combs[2,i]]
      datalist[[i]] <- c(x, y, 
                       round(RMSE(my.dat[,x],my.dat[, y]), 3),
                       sum(!is.na(my.dat[,x] & my.dat[,y])))
    }
  
    rmse_data <- as.data.frame(do.call(rbind, datalist))
    names(rmse_data) <- c("Lib1", "Lib2", "rmse", "nsites")
    my.rmse[[tmp.name]] <- rmse_data
  }
  
  # 6 plotting
  #Randomly sample 10k CGs to plot the smooth scatter: otherwise it takes too long to plot
  {
    print("plotting")
    df <- as.data.frame(my.dat[sample(nrow(my.dat), 10000), ])

    p <- ggpairs(df, lower="blank", upper = "blank") 
    seq <- 1:ncol(df)
  
    for (x in seq)
      for (y in seq) 
        if (y>x){
          p <- putPlot(p, ggplot(df, aes_string(x=names(df)[x],y=names(df)[y])) + 
                         stat_density2d(aes(fill = ..density..^0.55), geom = "tile", contour = FALSE, n = 200) +
                         theme(axis.text.x = element_text(angle = 90))+
                         scale_fill_viridis_c(option = "plasma"), y,x)
        } 
        

   for (x in seq) 
      for (y in seq) 
        if (x>y) {
          rmse.test <- as.character(rmse_data$rmse[rmse_data$Lib1 %in% names(df)[x] & rmse_data$Lib2 %in% names(df)[y] |
                   rmse_data$Lib2 %in% names(df)[x] & rmse_data$Lib1 %in% names(df)[y]])
          n.sites <- as.character(rmse_data$nsites[rmse_data$Lib1 %in% names(df)[x] & rmse_data$Lib2 %in% names(df)[y] |
                  rmse_data$Lib2 %in% names(df)[x] & rmse_data$Lib1 %in% names(df)[y]])
          tmp.my.cor.p <- my.cor.p[[tmp.name]]
          tmp.my.cor.s <- my.cor.s[[tmp.name]]
          tmp.p <- tmp.my.cor.p$value[tmp.my.cor.p$Var1 %in% names(df)[x] & tmp.my.cor.p$Var2 %in% names(df)[y]]
          tmp.s <- tmp.my.cor.s$value[tmp.my.cor.s$Var1 %in% names(df)[x] & tmp.my.cor.s$Var2 %in% names(df)[y]]
     
        
         p <- putPlot(p, ggplot(df, aes_string(x=names(df)[x],y=names(df)[y])) + 
                     theme(panel.background = element_blank(), 
                           axis.title = element_blank(),
                           axis.text = element_blank(), 
                           axis.ticks = element_blank()) +
                     annotate("text", x=1, y=2,size =3, 
                              label= sprintf("r= %s \np= %s \nrmse= %s \n%s CpGs",
                                             tmp.p, tmp.s, rmse.test, n.sites)), y,x)
         }

    

    png(filename= sprintf("%s_AssayComparison.png", tmp.name), 
        width=1300, height=850, res=200, bg = "transparent")
     p
    dev.off()
  }

}


save(my.stats,my.cpgs, my.cor.p, my.cor.s, my.rmse, large.df, file = "Assay_compare_data.RData")


