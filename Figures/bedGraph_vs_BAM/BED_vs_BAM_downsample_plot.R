#Plotting the downsampling comparison
#Supplemental Fig S5
library(GGally)
library(data.table)


# FUNCTIONS #
RMSE = function(m, o){
  sqrt(mean((m - o)^2, na.rm=T))
}


# Load file with Beta-values for all CpG sites >=5x coverage from bedGraph or BAM downsampling or original BAM file for HG006 EMSeq lab 1
my.dat <- fread('beta_HG006_LAB01_5x.txt.gz')
my.dat <- my.dat[,c(7,5,8,6,9)] # select only bedGraphs downsamples with https://github.com/nebiolabs/methylation_tools/downsample_methylKit.py method
names(my.dat) <- c( "BAM_10x", "bedGraph_10x", "BAM_20x", "BedGraph_20x", "BAM_44x") #simplify names


# Calculate genome-wide correlations, RMSE and Pair-wise overlapping sites
count.cpgs.covered <- colSums(!is.na(my.dat)) #COunt CPG sites 
my.cor.p <- round(cor(my.dat, use = "pairwise.complete.obs", method="pearson"), 3)
my.cor.s <- round(cor(my.dat, use = "pairwise.complete.obs", method="spearman"), 3)

cor.p <- reshape2::melt(my.cor.p)
cor.s <- reshape2::melt(my.cor.s)


# Calc stats for each library
stats.beta <- sapply(my.dat, function(x) summary(x))
stats.beta <- rbind(stats.beta, count.cpgs.covered)
write.table(stats.beta, file="EMSeq_HG006_DownSample_BetaValue_statistics.txt", sep="\t", quote=F, row.names=T)

# Calc RMSE
combs <- combn(1:ncol(my.dat), 2) 
datalist = list()
for (i in 1:ncol(combs)) {
  x <- names(my.dat)[combs[1,i]]
  y <- names(my.dat)[combs[2,i]]
  datalist[[i]] <- c(x, y, 
                     round(RMSE(my.dat[,get(x)],my.dat[, get(y)]), 3),
                     sum(!is.na(my.dat[,get(x)] & my.dat[,get(y)])))
}

# collect & organize RMSE, Pearsons and Spearmans correlation values for plots
rmse_data = as.data.frame(do.call(rbind, datalist))
names(rmse_data) <- c("Lib1", "Lib2", "rmse", "nsites")
rmse_data$name_merge <- paste(rmse_data$Lib1, rmse_data$Lib2, sep="_")
cor.p$name_merge <- paste(cor.p$Var1, cor.p$Var2, sep="_")
cor.s$name_merge <- paste(cor.s$Var1, cor.s$Var2, sep="_")
test <- merge.data.frame(rmse_data, cor.p, by= "name_merge")
my.merged.df <- merge.data.frame(test, cor.s, by="name_merge")
names(my.merged.df)[8] <- "pearson_rho"
names(my.merged.df)[11] <- "spearman_p"
my.merged.df$nsites <- round(as.numeric(as.character(my.merged.df$nsites))/1000000, 2)



# Plot comparison of the orginal bedGraph (44x BAM) to the downsampled data
m.df <- reshape2::melt(my.merged.df, id.vars = c("name_merge", "nsites"), measure.vars= c("rmse", "pearson_rho", "spearman_p"))
m.df$value <- as.numeric(m.df$value)
v44x <-  m.df[grepl("44x", m.df$name_merge),]
v44x$order <- rep(c(1,3,2,4), 3)
v44x$order <- rep(c("10x BAM",  "20x BAM", "10x bedGraph","20x bedGraph"), 3)
ggplot(data = v44x, aes(order, value)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color="steelblue") + 
  labs(title = "",
       subtitle = "Metrics by downsampling method in comparison to original 44x BAM",
       y = "metric", x = "") + 
  facet_wrap(~ variable,scales = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

 
#Randomly sample 10k CGs to plot the histogram and smooth scatter: otherwise it takes too long
df <- as.data.frame(my.dat[sample(nrow(my.dat), 10000), ])

# MAKE A SCATTER PLOT MATRIX
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
p



###############
## Plot count coverage distributions in bedgraphs only for chr 1 (to save computational time)
dat.dir <- 'chr1_bedgraphs'
my.files <- list.files(dat.dir, recursive=F, full.names= T, pattern= "chr1.gz")

my.stats <- as.list(NULL)
hist.dat <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("name", "CG_coverage", "assay"))

# load & process counts for chr 1
for (f in my.files) {
  tmp <- fread(f)
  #fix some names
  tmp.name <- sub("_bedGraph_chr1.gz", "", basename(f))
  tmp.name <- sub("_na", "", tmp.name)
  tmp.assay <- unlist(strsplit(tmp.name, "_"))[1]
  tmp.name.plot <- paste(unlist(strsplit(tmp.name, "_"))[-1], collapse ="_")
  #get the counts
  counts <- tmp[,get("V5")] + tmp[,get("V6")]
  #basic stats
  my.stats[[tmp.name]] <- c(summary(counts), quantile(counts, .99), quantile(counts, .95))
  #data for making histograms
  hist.dat <- rbind(hist.dat, data.frame(name= rep(tmp.name.plot, length(counts)), 
                                         CG_coverage=counts, 
                                         assay=rep(tmp.assay, length(counts))))
}

df1 <-as.data.frame(my.stats)

# Plot a coverage histogram of
emseq <- hist.dat[hist.dat$assay %in% "EMSeq" ,]
emseq <- emseq[!grepl("TotC", emseq$name) ,] #remove totC, ie only keep the https://github.com/nebiolabs/methylation_tools/downsample_methylKit.py downsampling method for bedGraphs
emseq$new.name <- as.factor(as.character(emseq$name))
emseq$new.name <- factor(emseq$new.name,levels(emseq$new.name)[c(3,6,1,4,7,2,5)])
ggplot(emseq, aes(x=CG_coverage, fill=new.name)) +
  geom_histogram(alpha=0.4, binwidth = 1, position='identity') +
  xlim(0,75) +
  theme(legend.title=element_blank())

