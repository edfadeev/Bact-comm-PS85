#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#load libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("ggpmisc"); packageVersion("ggpmisc")
library("splines"); packageVersion("splines")
library("vegan"); packageVersion("vegan")
library("cowplot"); packageVersion("cowplot")

#set plots theme
theme_set(theme_classic())

#load colour palettes
source('./Scripts/color_palettes.R')

#####################################
#Figure 5: Bacterial community characteristics across the Fram Strait
####################################
BAC_pruned.vst <- readRDS("./Data/BAC_pruned_vst.rds")

#generate data set with all bacterial community characteristics
BAC_comm.char<- data.frame(Expedition = sample_data(BAC_pruned.vst)$Expedition,
                           Longitude = sample_data(BAC_pruned.vst)$Longitude..degrees_east.,
                           Type = sample_data(BAC_pruned.vst)$Type,
                           Fraction = sample_data(BAC_pruned.vst)$Fraction,
                           Region = sample_data(BAC_pruned.vst)$Region,
                           Cell.conc = log10(sample_data(BAC_pruned.vst)$Conc_cells.mL.),
                           HNA.LNA = sample_data(BAC_pruned.vst)$FCM_HNA/sample_data(BAC_pruned.vst)$FCM_LNA,
                           Bac.prod. = sample_data(BAC_pruned.vst)$leucine,
                           spec.bac.prod = (sample_data(BAC_pruned.vst)$leucine)/(1000*sample_data(BAC_pruned.vst)$Conc_cells.mL.))

#cell concentration plot
BAC_cell.conc.p <-   ggplot(BAC_comm.char[!is.na(BAC_comm.char$Cell.conc),],
                            aes(x = Longitude, y = Cell.conc,  group = Expedition, colour = Region, shape = Type)) + #add 'group = Type' for smooth line
  geom_point(size=4) + 
  xlab("Longitude [°East]")+
  ylab("Cell density [log10(cells/ml)]")+
  theme(legend.position="none")+
  scale_shape_manual(values=c(3, 4, 7))+
  geom_smooth(aes(linetype = Expedition), method = "gam", show.legend = TRUE,formula =  y ~ ns(x,4), se=FALSE, colour = "black")+
  scale_colour_manual(values=reg_colours)+
  stat_poly_eq(formula = y ~ ns(x,4), aes(label = ..rr.label..), parse = TRUE)+ 
  xlim(-11,6)

#bacterial productivity plot
BAC_prod.p <-   ggplot(BAC_comm.char[!is.na(BAC_comm.char$Bac.prod.),],
                       aes(x = Longitude, y = Bac.prod., group = Expedition, colour = Region, shape = Type)) + #add 'group = Type' for smooth line
  ylab("Bacterial production [pmol leucine L-1 h-1]")+
  geom_point(size=4) + 
  theme(legend.position="none")+
  xlab("Longitude [°East]")+
  scale_shape_manual(values=c(3, 4, 7))+
  geom_smooth(aes(linetype = Expedition), method = "gam", show.legend = TRUE,formula =  y ~ ns(x,4), se=FALSE, colour = "black")+
  scale_colour_manual(values=reg_colours)+
  stat_poly_eq(formula = y ~ ns(x,4),aes(label = ..rr.label..),parse = TRUE)+
  xlim(-11,6)

#####################################
#PCoA plot
#####################################
BAC_comm.ord <- ordinate(BAC_pruned.vst, method = "MDS", distance = "euclidean")
evals <- BAC_comm.ord$values$Eigenvalues
BAC_comm.ord.p <- plot_ordination(BAC_pruned.vst, BAC_comm.ord, color = "Region", shape = "Fraction") +
  labs(fill = "Water layer", shape = "Fraction") +
  scale_colour_manual(values=c("EGC"="blue", "WSC"="red")) +
  geom_point(aes(fill = Region, shape = Fraction), size = 4) +
  geom_point(fill = "grey90", size = 1.5)+
  geom_text(aes(label = StationName), nudge_y= -5,  size=3)+
  stat_ellipse(geom = "polygon", alpha = 0.1)+
  theme(legend.position="none")+ 
  scale_x_reverse()+
  scale_y_reverse()

#combined plot
png('./Figures/Figure-5-bacterial-overview.png',width = 30, height = 30, res = 300,units = "cm")
plot_grid(BAC_cell.conc.p, BAC_prod.p, BAC_comm.ord.p, labels = c("A", "B", "C"), ncol = 2, align = "hv")
dev.off()

#check whether the regions definition is significant
BAC_metadata <- as(sample_data(BAC_pruned.vst), "data.frame")
BAC_comm.dist <- phyloseq::distance(BAC_pruned.vst, "euclidean")
BAC_comm.adonis <- adonis(BAC_comm.dist ~ Region + Fraction + Type , BAC_metadata)
BAC_comm.adonis

#####################################
#PCoA Eukaryotes
####################################
EUK_pruned.vst <- readRDS("./Data/EUK_pruned_vst.rds")

EUK_comm.ord <- ordinate(EUK_pruned.vst, method = "MDS", distance = "euclidean")
evals_EUK <- EUK_comm.ord$values$Eigenvalues
EUK_comm.ord.p <- plot_ordination(EUK_pruned.vst, EUK_comm.ord, color = "Region") +
  labs(fill = "Water layer", shape = "Fraction") +
  stat_ellipse(geom = "polygon", alpha = 0.1)+
  scale_colour_manual(values=c("EGC"="blue", "WSC"="red")) +
  geom_point(aes(fill = Region, group = Region, shape = Type), size = 6) +
  geom_text(aes(label = StationName), colour = "black",nudge_y= -5,  size=6)+
  theme(legend.position="none")+ 
  scale_x_reverse()

#check whether the regions definition is significant
EUK_metadata <- as(sample_data(EUK_pruned.vst), "data.frame")
EUK_comm.dist <- phyloseq::distance(EUK_pruned.vst, "euclidean")
EUK_comm.adonis <- adonis(EUK_comm.dist ~ Region + Type , EUK_metadata)
EUK_comm.adonis

#####################################
#remove all temporary datasets 
####################################
detach("package:phyloseq", unload=TRUE)
detach("package:ggpmisc", unload=TRUE)
detach("package:ggplot2", unload=TRUE)
detach("package:DESeq2", unload=TRUE)
detach("package:cowplot", unload=TRUE)
detach("package:vegan", unload=TRUE)

rm(list=ls(all=TRUE))

sessionInfo()
