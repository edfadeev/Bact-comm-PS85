#####################################
#working directory
#####################################

setwd("C:/Users/Eddie/Google Drive/PhD/Chapter_1_Spacial_comparison/16S_analysis/surface_water_analysis_in_R/")

#####################################
#Load libraries
#####################################

library("phyloseq")
library("ggplot2")
library("DESeq2")
library("vegan")
source("../OTU_analysis/multiplot.R")

#load image
load(file="./phyloseq.RData")

#######Select the proper dataframe according to data_import.R
# X1 =  subset of the dataset without any manipulations
# X1_red = dataset after prevalence filter
# X1_sub = subset of entries with known taxonomic 
# X1_rare = subset of rarified samples, according to the smallest library
# X1_Hvar = subset of OTU with high variation
# X1_log = log transformed counts

####based on http://joey711.github.io/phyloseq-extensions/DESeq2.html

#####################################
#Subset the samples into different two groups
#####################################

#group by cruise
X1_rare.sub <- subset_samples(X1_rare, Expedition=="PS85" | Expedition=="PS85")

#group by water mass
X1_rare.sub <- subset_samples(X1_rare, Water.mass=="PSW" | Water.mass=="AW")

#group by sample type
X1_rare.sub <- subset_samples(X1_rare, Type== "SRF" | Type== "DCM")

#group by station name
X1_rare.sub <- subset_samples(X1_rare, StationName=="HG9" | StationName=="HG1")

#nested subset by sample type
X1_rare.sub <- subset_samples(subset_samples(X1_rare, Water.mass=="PSW" | Water.mass=="AW"), Expedition=="PS99.2")


#####################################
#Analyze phyloseq dataset in DESeq2
#####################################

#Define grouping using 'design' argument
X1_rare.DEseq <- phyloseq_to_deseq2(X1_rare.sub, design= ~Water.mass)

#define test
X1_rare.DEseq <- DESeq(X1_rare.DEseq, test="Wald", fitType="parametric")

#####################################
#Investigate test results
#####################################
X1_rare.DEseq.res <- results(X1_rare.DEseq, cooksCutoff = FALSE)

#define alpha
alpha <- 0.001

X1_rare.DEseq.sig <- X1_rare.DEseq.res[which(X1_rare.DEseq.res$padj < alpha), ]
X1_rare.DEseq.sig <- cbind(as(X1_rare.DEseq.sig, "data.frame"), as(tax_table(X1_rare.sub)[rownames(X1_rare.DEseq.sig), ], "matrix"))
head(X1_rare.DEseq.sig)

#####################################
#Visualize the results
#####################################
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x <- tapply(X1_rare.DEseq.sig$log2FoldChange, X1_rare.DEseq.sig$Class, function(x) max(x))
x <- sort(x, TRUE)
X1_rare.DEseq.sig$Class <- factor(as.character(X1_rare.DEseq.sig$Class), levels=names(x))

# Order order
x <- tapply(X1_rare.DEseq.sig$log2FoldChange, X1_rare.DEseq.sig$Order, function(x) max(x))
x <- sort(x, TRUE)
X1_rare.DEseq.sig$Order <- factor(as.character(X1_rare.DEseq.sig$Order), levels=names(x))

windows()
ggplot(X1_rare.DEseq.sig, aes(x=Order, y=log2FoldChange, color=Class)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=12),
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5),
              panel.background = element_rect(fill = 'white', colour = 'white'),
              panel.margin=unit(.05, "lines"),
              strip.text = element_text(size = 15),
              panel.border = element_rect(color = "black", fill = NA, size = 1), 
              strip.background = element_rect(color = "black", size = 1),
              legend.text= element_text(size = 12),
              legend.title= element_text(size = 14))+
  geom_hline(aes(yintercept=0), linetype="dashed")
