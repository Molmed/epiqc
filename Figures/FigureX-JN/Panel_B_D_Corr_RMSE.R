# Code for Panel B of main figure
# Jessica Nordlund
# June 2020

library(ggplot2) 
library(ggbeeswarm)
library(corrplot)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)

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



setwd("/Users/jessicanordlund/Downloads/EPIQC/for_manus")
load("Assay_compare_data.RData")

#####################################################
# Plot Correlation data
#####################################################

cor.p.final <- as.data.frame(lapply(my.cor.p, "[", , "value"))
cor.p.final$assay <- as.character(lapply(my.cor.p, "[", , "Var2")[[1]])

row.has.one <- apply(cor.p.final, 1, function(x){any(x == "1.000")})
row.has.one[row.has.one %in% NA] <- FALSE

cor.p.final <- cor.p.final[!row.has.one ,]
cor.p.final$value.type <- rep("rho", dim(cor.p.final)[1])

cor.p.final <- reshape2::melt(cor.p.final, id.vars = c("assay", "value.type"))
cor.p.final$value <- as.numeric(as.character(cor.p.final$value))
names(cor.p.final)[4] <- "pearsons_r"
names(cor.p.final)[3] <- "genome"

plot_cor = ggplot(cor.p.final, aes(x=assay, y=pearsons_r, color=assay)) +
  geom_boxplot(show.legend=F, aes(fill=assay), alpha=0.3) +
  geom_beeswarm(cex=2, size=2, show.legend=T, aes(shape=genome)) +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=0, size=12, hjust=0.5, vjust=0.5),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("EMSeq"="#4daf4a",
                                "MethylSeq"="#377eb8", 
                                "SPLAT"="#984ea3",
                                "TrueMethyl" = "#a65628", 
                                "TruSeq"= "#e41a1c")) + 
  scale_fill_manual(values  = c("EMSeq"="#4daf4a",
                                "MethylSeq"="#377eb8", 
                                "SPLAT"="#984ea3",
                                "TrueMethyl" = "#a65628", 
                                "TruSeq"= "#e41a1c")) +
  scale_shape_manual(values=c(0, 1, 2, 3, 10, 15, 16, 17 ))

png(filename= "AssayCompare_Pearsons_box_10x.png", width= 1300, heigh=900, res=200,bg = "transparent") 
plot_cor
dev.off()



#################################################
# PLOT RMSE data
#################################################
rmse.final <- as.data.frame(lapply(my.rmse, "[", , "rmse"))
rmse.final <- rbind(rmse.final, rmse.final)
rmse.final$assay <- c(as.character(lapply(my.rmse, "[", , "Lib1")[[1]]), as.character(lapply(my.rmse, "[", , "Lib2")[[1]]))
rmse.final$value.type <- rep("rmse", dim(rmse.final)[1])

rmse.final <- reshape2::melt(rmse.final, id.vars = c("assay", "value.type"))
rmse.final$value <- as.numeric(as.character(rmse.final$value))
rmse.final$value[rmse.final$value %in% "NaN"] <- NA
names(rmse.final)[4] <- "rmse"
names(rmse.final)[3] <- "genome"

plot_cor = ggplot(rmse.final, aes(x=assay, y=rmse, color=assay)) +
  geom_boxplot(show.legend=F, aes(fill=assay), alpha=0.3) +
  geom_beeswarm(cex=2, size=2, show.legend=T, aes(shape=genome)) +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=0, size=12, hjust=0.5, vjust=0.5),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("EMSeq"="#4daf4a",
                                "MethylSeq"="#377eb8", 
                                "SPLAT"="#984ea3",
                                "TrueMethyl" = "#a65628", 
                                "TruSeq"= "#e41a1c")) + 
  scale_fill_manual(values  = c("EMSeq"="#4daf4a",
                                "MethylSeq"="#377eb8", 
                                "SPLAT"="#984ea3",
                                "TrueMethyl" = "#a65628", 
                                "TruSeq"= "#e41a1c")) +
  scale_shape_manual(values=c(0, 1, 2, 3, 10, 15, 16, 17 ))

png(filename= "AssayCompare_RMSE_box_10x.png", width= 1300, heigh=900, res=200,bg = "transparent") 
plot_cor
dev.off()




#################################################
# PLOT PAIRWISE CORRELATIONS
#################################################
c.df <- round(cor(large.df, use = "pairwise.complete.obs", method="pearson"),3)

make.labels <- matrix(unlist(strsplit( colnames(c.df), "_")), ncol=2, byrow=TRUE)
top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
                                                    labels = c("group1", "group2", "group3"), 
                                                    labels_gp = gpar(col = "white", fontsize = 10)))

ha = HeatmapAnnotation(genome = make.labels[,1],
                       assay= make.labels[,2],
                       col = list(genome= my.genome.cols, 
                                  assay=my.assay.cols))
hb = rowAnnotation(genome = make.labels[,1],
                       assay= make.labels[,2],
                       col = list(genome= my.genome.cols, 
                                  assay=my.assay.cols))

hc = HeatmapAnnotation(foo = anno_empty(border = TRUE))

Heatmap(c.df, name = "Pearson's r", top_annotation = ha, 
        left_annotation = hb, show_row_names = FALSE, show_column_names = FALSE)

