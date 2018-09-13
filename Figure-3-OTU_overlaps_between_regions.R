#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#load libraries
library(phyloseq); packageVersion("phyloseq")
library(UpSetR); packageVersion("UpSetR")
library(ggplot2); packageVersion("ggplot2")
library(gridExtra); packageVersion("gridExtra")
library(VennDiagram); packageVersion("VennDiagram")

source("./Scripts/pres_abs_matrix.R")
source("./Scripts/color_palettes.R")

#####################################
#Diagram of bacterial OTU overlap between the different groups 
#####################################
BAC_pruned <- readRDS("./Data/BAC_pruned.rds")

#FL
BAC_FL <- subset_samples(BAC_pruned, Fraction == "0.22")
BAC_FL <- prune_taxa(taxa_sums(BAC_FL)>0,BAC_FL)
BAC_FL.WSC <- subset_samples(BAC_FL, Region == "WSC")
BAC_FL.WSC <- prune_taxa(taxa_sums(BAC_FL.WSC)>0,BAC_FL.WSC)
BAC_FL.EGC <- subset_samples(BAC_FL, Region == "EGC")
BAC_FL.EGC <- prune_taxa(taxa_sums(BAC_FL.EGC)>0,BAC_FL.EGC)
#PA
BAC_PA <- subset_samples(BAC_pruned, Fraction == "3")
BAC_PA <- prune_taxa(taxa_sums(BAC_PA)>0,BAC_PA)
BAC_PA.WSC <- subset_samples(BAC_PA, Region == "WSC")
BAC_PA.WSC <- prune_taxa(taxa_sums(BAC_PA.WSC)>0,BAC_PA.WSC)
BAC_PA.EGC <- subset_samples(BAC_PA, Region == "EGC")
BAC_PA.EGC <- prune_taxa(taxa_sums(BAC_PA.EGC)>0,BAC_PA.EGC)

#make a list
y <- list()
y[["WSC-FL"]] <- as.character(row.names(otu_table(BAC_FL.WSC)))
y[["WSC-PA"]] <- as.character(row.names(otu_table(BAC_PA.WSC)))
y[["EGC-FL"]] <- as.character(row.names(otu_table(BAC_FL.EGC)))
y[["EGC-PA"]] <- as.character(row.names(otu_table(BAC_PA.EGC)))

#generate overlap matrix
otu_overlaps <- pres_abs_matrix(y)    
otu_overlaps$OTU <- rownames(otu_overlaps)

#####################################
#Explore overlaping OTU
#####################################
taxonomy <- as.data.frame(tax_table(BAC_pruned))
taxonomy$OTU <- rownames(taxonomy)

otu_overlaps_merged <- full_join(taxonomy,otu_overlaps, by = c("OTU"))

#####################################
#Relative abundance of overlaping OTU
#####################################
#add abundance in each fraction
#transform data
BAC_pruned.ra <- transform_sample_counts(BAC_pruned, function(x) x / sum(x))

#calculate mean abundance for each OTU
BAC_pruned.ra.long <- psmelt(BAC_pruned.ra)
BAC_pruned.ra.long.agg <- aggregate(Abundance~OTU+Fraction, BAC_pruned.ra.long, FUN = mean)
BAC_pruned.ra.long.agg$Abundance <- BAC_pruned.ra.long.agg$Abundance*100

otu_overlaps_merged <- full_join(otu_overlaps_merged, BAC_pruned.ra.long.agg[BAC_pruned.ra.long.agg$Fraction=="0.22",], by = "OTU")

otu_overlaps_merged <- full_join(otu_overlaps_merged, BAC_pruned.ra.long.agg[BAC_pruned.ra.long.agg$Fraction=="3",], by = "OTU")

colnames(otu_overlaps_merged)[13] <- "Abundance.FL"
colnames(otu_overlaps_merged)[15] <- "Abundance.PA"

#set metadata
sets <- names(y)
metadata <- as.data.frame(cbind(sets, rep(c("FL", "PA"),2),
                                c(rep("WSC",2),rep("EGC",2))))
names(metadata) <- c("sets", "Fraction","Region")

#plot
png('./Figures/Figure-3-OTU-overlap.png',width = 30, height = 30, units = "cm")
upset(otu_overlaps_merged, number.angles = 30,
      sets = as.vector(metadata$sets),
      keep.order = TRUE, 
      mainbar.y.label = "No. of overlaping OTU",
      order.by = "freq", #empty.intersections = "on",
      mainbar.y.max = 1500,
      #group.by = "degree",
      boxplot.summary = c("Abundance.FL", "Abundance.PA"),
      queries = list(list(query = intersects, 
                          params = list("WSC-FL","EGC-FL","WSC-PA","EGC-PA"), 
                          color = "yellow"),
                     list(query = intersects, 
                          params = list("WSC-FL","WSC-PA"), 
                          color = "red"),
                      list(query = intersects, 
                           params = list("EGC-FL","EGC-PA"), 
                           color = "blue")),
      set.metadata = list(data = metadata, plots = list(list(type = "matrix_rows", 
                                                             column = "Region", colors = c(EGC = "blue", WSC = "red"), 
                                                             alpha = 0.5))))

dev.off()
      
#plot relative abundance     
BAC_pruned.ra.long.shared <- BAC_pruned.ra.long[BAC_pruned.ra.long$OTU %in% rownames(pres_abs_matrix(y))[rowSums(pres_abs_matrix(y))==4],]

levels(BAC_pruned.ra.long.shared$StationName) <- c("10W","8.5W","7W","EG1","EG3","EG4",
                                        "1W","1E","HG9","N4","HG4","HG1")

#aggregate by taxonomy
BAC_pruned.ra.long.shared.agg <- aggregate(Abundance~StationName+Type+Fraction+Class, BAC_pruned.ra.long.shared, FUN= sum)

#remove beloew 1% ra
BAC_pruned.ra.long.shared.agg$Abundance <- BAC_pruned.ra.long.shared.agg$Abundance*100
BAC_pruned.ra.long.shared.agg <- BAC_pruned.ra.long.shared.agg[BAC_pruned.ra.long.shared.agg$Abundance>0.5,]

  
#define number of colours
n_shared <- length(levels(droplevels(BAC_pruned.ra.long.shared.agg$Class)))
cols_shared <- sample(tol21rainbow, n_shared, replace= TRUE)

#plot
BAC_shared.otu.plot <- ggplot(BAC_pruned.ra.long.shared.agg, aes(x = StationName, y = Abundance, fill = Class)) + 
  facet_grid(Type~Fraction, space= "fixed") +
  geom_col()+
  scale_fill_manual(values = tol21rainbow) +
  theme(legend.position = "bottom")+ 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance of shared OTU (Class > 0.5%) \n")

ggsave("./Figures/Supp-Figure-2-RA_shared_OTU.png", BAC_shared.otu.plot,
       dpi = 300, device= "png",width = 30, height = 30, units = "cm")


#####################################
#Venn diagram of Eukaryotic OTU overlap between the different regions 
#####################################
EUK_pruned <- readRDS("./Data/EUK_pruned.rds")

EUK_WSC <- subset_samples(EUK_pruned, Region == "WSC")
EUK_WSC <- prune_taxa(taxa_sums(EUK_WSC)>0, EUK_WSC)

EUK_EGC <- subset_samples(EUK_pruned, Region == "EGC")
EUK_EGC <- prune_taxa(taxa_sums(EUK_EGC)>0, EUK_EGC)


EUK_over <- list()
EUK_over$WSC <- taxa_names(EUK_WSC)
EUK_over$EGC <- taxa_names(EUK_EGC)

EUK_over.venn <- venn.diagram(EUK_over,
  lwd = 1, fill = c("red", "yellow"),  alpha = c(0.5, 0.5), cex = 2,
  #cat.fontface = 4,
  lty =2, filename = NULL, scaled = TRUE,
  inverted = FALSE, print.mode = c("raw","percent"))

grid.newpage()
grid.draw(EUK_over.venn)

EUK_over.calc <- calculate.overlap(EUK_over)

#transform data
EUK_pruned.ra <- transform_sample_counts(EUK_pruned, function(x) x / sum(x))

#calculate mean abundance for each OTU
EUK_pruned.ra.long <- psmelt(EUK_pruned.ra)

#plot relative abundance     
EUK_pruned.ra.long.shared <- EUK_pruned.ra.long[EUK_pruned.ra.long$OTU %in% EUK_over.calc$a3,]

levels(EUK_pruned.ra.long.shared$StationName) <- c("10W","8.5W","7W","EG1","EG3","EG4",
                                                   "1W","1E","HG9","N4","HG4","HG1")

#define number of colours
n_shared <- length(levels(EUK_pruned.ra.long.shared$X3))
cols_shared <- sample(tol21rainbow, n_shared, replace= FALSE)

#plot
EUK_pruned.ra.long.shared.p <- ggplot(EUK_pruned.ra.long.shared, aes(x = StationName, y = Abundance, fill = X3)) + 
  facet_grid(Type~., space= "fixed") +
  geom_col()+
  scale_fill_manual(values = tol21rainbow) +
  theme(legend.position = "bottom")+ 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance of shared OTU (%) \n")

ggsave("./Figures/Supp-Figure-3-Overlaps-euk.png", EUK_pruned.ra.long.shared.p,
       dpi = 300, device= "png",width = 30, height = 30, units = "cm")

#####################################
#remove all temporary datasets 
####################################
detach("package:phyloseq", unload=TRUE)
detach("package:UpSetR", unload=TRUE)
detach("package:ggplot2", unload=TRUE)
detach("package:gridExtra", unload=TRUE)
detach("package:VennDiagram", unload=TRUE)

rm(list=ls(all=TRUE))

sessionInfo()
