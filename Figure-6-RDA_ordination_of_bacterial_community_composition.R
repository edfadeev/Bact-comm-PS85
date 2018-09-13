#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#load libraries
library("ggplot2"); packageVersion("ggplot2")
library("ggpmisc"); packageVersion("ggpmisc")
library("phyloseq"); packageVersion("phyloseq")
library("vegan"); packageVersion("vegan")
library("cowplot"); packageVersion("cowplot")
library("dplyr"); packageVersion("dplyr")

#set plots theme
theme_set(theme_classic())

#load colour palettes
source('./Scripts/color_palettes.R')

#####################################
#Figure 6: RDA ordination of the fitted model of bacterial composition data constrained by environmental variables
####################################
BAC_pruned.vst <- readRDS("./Data/BAC_pruned_vst.rds")

#Free-living fraction
#subset by fraction and remove NAs in metadata
BAC_FL.no.na <- BAC_pruned.vst %>% 
  subset_samples(
    !is.na(dNO3) & 
      !is.na(dSiO3) &
      !is.na(dPO4) &
      !is.na(ChlA)&
      Fraction =="0.22")
##remove unobserved OTU
BAC_FL.no.na <- prune_taxa(taxa_sums(BAC_FL.no.na)>0,BAC_FL.no.na)

#extract and scale the env. parameters
BAC_FL.env <- data.frame(sample_data(BAC_FL.no.na))[c("Temperature", "Salinity", "dNO3", "dSiO3", "dPO4", "ChlA")]  
BAC_FL.env <- as.data.frame(scale(BAC_FL.env,center = FALSE, scale = TRUE))

#extract OTU tables from Phyloseq object
BAC_FL.otu <- t(otu_table(BAC_FL.no.na))

#Particles-associated fraction
BAC_PA.no.na <- BAC_pruned.vst %>% 
  subset_samples(
    !is.na(dNO3) & 
      !is.na(dSiO3) &
      !is.na(dPO4) &
      !is.na(ChlA)&
      Fraction =="3")
##remove unobserved OTU
BAC_PA.no.na <- prune_taxa(taxa_sums(BAC_PA.no.na)>0,BAC_PA.no.na)

#extract and scale the env. parameters
BAC_PA.env <- data.frame(sample_data(BAC_PA.no.na))[c("Temperature", "Salinity", "dNO3", "dSiO3", "dPO4", "ChlA")]  
BAC_PA.env <- as.data.frame(scale(BAC_PA.env,center = FALSE, scale = TRUE))

#extract OTU tables from Phyloseq object
BAC_PA.otu <- t(otu_table(BAC_PA.no.na))


#RDA analysis
BAC_FL.rda.all <- rda (BAC_FL.otu ~ ., data = BAC_FL.env) # model including all variables 
BAC_PA.rda.all <- rda (BAC_PA.otu ~ ., data = BAC_PA.env) # model including all variables 

#generate an RDA plot 
#FL
BAC_FL.rda.scores <- vegan::scores(BAC_FL.rda.all,display=c("sp","wa","lc","bp","cn"))
BAC_FL.rda.sites <- data.frame(BAC_FL.rda.scores$sites)
BAC_FL.rda.sites$Sample.ID <- as.character(rownames(BAC_FL.rda.sites))
sample_data(BAC_FL.no.na)$Sample.ID <- as.character(rownames(sample_data(BAC_FL.no.na)))
sample_data(BAC_FL.no.na)$Type <- as.character(sample_data(BAC_FL.no.na)$Type )
BAC_FL.rda.sites <- BAC_FL.rda.sites %>%
  left_join(sample_data(BAC_FL.no.na))

#Draw biplots
BAC_FL.rda.arrows<- BAC_FL.rda.scores$biplot*5
colnames(BAC_FL.rda.arrows)<-c("x","y")
BAC_FL.rda.arrows <- as.data.frame(BAC_FL.rda.arrows)
BAC_FL.rda.evals <- 100 * (BAC_FL.rda.all$CCA$eig / sum(BAC_FL.rda.all$CCA$eig))

#Plot 
BAC_FL.rda.plot <- ggplot() +
  geom_point(data = BAC_FL.rda.sites, aes(x = RDA1, y = RDA2, fill = Region), 
             shape =21, size = 4) +
  geom_text(data = BAC_FL.rda.sites,aes(x = RDA1, y = RDA2,label = StationName), 
            nudge_y= -0.3,size=3)+
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_FL.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_FL.rda.evals[2], 2))) +
  scale_fill_manual(values = reg_colours) +
  #scale_x_reverse()+ 
  theme(legend.position = "none")+
  geom_segment(data=BAC_FL.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(BAC_FL.rda.arrows*1.2),
            aes(x, y, label = rownames(BAC_FL.rda.arrows)),color="black",alpha=0.5)

#PA
BAC_PA.rda.scores <- vegan::scores(BAC_PA.rda.all,display=c("sp","wa","lc","bp","cn"))
BAC_PA.rda.sites <- data.frame(BAC_PA.rda.scores$sites)
BAC_PA.rda.sites$Sample.ID <- as.character(rownames(BAC_PA.rda.sites))
sample_data(BAC_PA.no.na)$Sample.ID <- as.character(rownames(sample_data(BAC_PA.no.na)))
sample_data(BAC_PA.no.na)$Type <- as.character(sample_data(BAC_PA.no.na)$Type )
BAC_PA.rda.sites <- BAC_PA.rda.sites %>%
  left_join(sample_data(BAC_PA.no.na))

#Draw biplots
BAC_PA.rda.arrows<- BAC_PA.rda.scores$biplot*5
colnames(BAC_PA.rda.arrows)<-c("x","y")
BAC_PA.rda.arrows <- as.data.frame(BAC_PA.rda.arrows)
BAC_PA.rda.evals <- 100 * (BAC_PA.rda.all$CCA$eig / sum(BAC_PA.rda.all$CCA$eig))

#Plot 
BAC_PA.rda.plot <- ggplot() +
  geom_point(data = BAC_PA.rda.sites, aes(x = RDA1, y = RDA2, fill = Region), 
             shape =21, size = 4) +
  geom_text(data = BAC_PA.rda.sites,aes(x = RDA1, y = RDA2,label = StationName), 
            nudge_y= -0.3,size=3)+
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_PA.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_PA.rda.evals[2], 2))) +
  scale_fill_manual(values = reg_colours) +
  #scale_x_reverse()+ 
  theme(legend.position = "none")+
  geom_segment(data=BAC_PA.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(BAC_PA.rda.arrows*1.2),
            aes(x, y, label = rownames(BAC_PA.rda.arrows)),color="black",alpha=0.5)


#combined figure
png('./Figures/Figure-6-RDA.png',width = 30, res = 300, height = 30, units = "cm")
plot_grid(BAC_FL.rda.plot, BAC_PA.rda.plot, labels = c("A", "B"), ncol = 2, align = "hv")
dev.off()

#####################################
#Forward selection of explanatory variables
#####################################
#based on the tutorial from David Zeleny Lab
#http://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel
#FL
BAC_FL.rda.0 <- rda (BAC_FL.otu ~ 1, data = BAC_FL.env) # model containing only species matrix and intercept
BAC_FL.rda.sel.os <- ordistep (BAC_FL.rda.0, scope = formula (BAC_FL.rda.all), direction = 'both') #stepwise selection

#PA
BAC_PA.rda.0 <- rda (BAC_PA.otu ~ 1, data = BAC_PA.env) # model containing only species matrix and intercept
BAC_PA.rda.sel.os <- ordistep (BAC_PA.rda.0, scope = formula (BAC_PA.rda.all), direction = 'both') #stepwise selection

#variance partitioning 
varpart(BAC_FL.otu, ~ Temperature, ~ ChlA, ~ Salinity, data = BAC_FL.env[,c("Temperature", "ChlA", "Salinity")])

varpart(BAC_PA.otu, ~ Salinity, ~ dNO3, ~ Temperature, ~ ChlA,data = BAC_PA.env[,c("Salinity", "Temperature", "ChlA","dNO3")])


#####################################
#Eukaryotic community RDA (supplementary)
####################################
EUK_pruned.vst <- readRDS("./Data/EUK_pruned_vst.rds")

#subset by fraction and remove NAs in metadata
EUK.no.na <- EUK_pruned.vst %>% 
  subset_samples(
    !is.na(dNO3) & 
      !is.na(dSiO3) &
      !is.na(dPO4) &
      !is.na(ChlA))
##remove unobserved OTU
EUK.no.na <- prune_taxa(taxa_sums(EUK.no.na)>0,EUK.no.na)

#extract and scale the env. parameters
EUK.env <- data.frame(sample_data(EUK.no.na))[c("Temperature", "Salinity", "dNO3", "dSiO3", "dPO4", "ChlA")]  
EUK.env <- as.data.frame(scale(EUK.env,center = FALSE, scale = TRUE))

#extract OTU tables from Phyloseq object
EUK.otu <- t(otu_table(EUK.no.na))


#RDA analysis
EUK.rda.all <- rda (EUK.otu ~ ., data = EUK.env) # model including all variables 

#generate an RDA plot 
EUK.rda.scores <- vegan::scores(EUK.rda.all,display=c("sp","wa","lc","bp","cn"))
EUK.rda.sites <- data.frame(EUK.rda.scores$sites)
EUK.rda.sites$Sample.ID <- as.character(rownames(EUK.rda.sites))
sample_data(EUK.no.na)$Sample.ID <- as.character(rownames(sample_data(EUK.no.na)))
sample_data(EUK.no.na)$Type <- as.character(sample_data(EUK.no.na)$Type )
EUK.rda.sites <- EUK.rda.sites %>%
  left_join(sample_data(EUK.no.na))

#Draw biplots
EUK.rda.arrows<- EUK.rda.scores$biplot*5
colnames(EUK.rda.arrows)<-c("x","y")
EUK.rda.arrows <- as.data.frame(EUK.rda.arrows)
EUK.rda.evals <- 100 * (EUK.rda.all$CCA$eig / sum(EUK.rda.all$CCA$eig))

#Plot 
EUK.rda.plot <- ggplot() +
  geom_point(data = EUK.rda.sites, aes(x = RDA1, y = RDA2, colour = Region, shape = Type), 
             size = 6) +
  geom_text(data = EUK.rda.sites,aes(x = RDA1, y = RDA2,label = StationName), 
            nudge_y= -0.3,size=6)+
  labs(x = sprintf("RDA1 [%s%%]", round(EUK.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(EUK.rda.evals[2], 2))) +
  scale_colour_manual(values = reg_colours) +
  scale_y_reverse()+ 
  theme(legend.position = "none")+
  geom_segment(data=EUK.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(EUK.rda.arrows*1.2),
            aes(x, y, label = rownames(EUK.rda.arrows)),color="black",alpha=0.5)


EUK.rda.0 <- rda (EUK.otu ~ 1, data = EUK.env) # model containing only species matrix and intercept
EUK.rda.sel.os <- ordistep (EUK.rda.0, scope = formula (EUK.rda.all), direction = 'both') #stepwise selection


#variance partitioning 
varpart(EUK.otu, ~ Temperature, ~ Salinity, ~ dSiO3, ~ dNO3, data = EUK.env[,c("Temperature", "Salinity", "dSiO3", "dNO3")])


#####################################
#remove all temporary datasets 
####################################
detach("package:phyloseq", unload=TRUE)
detach("package:ggpmisc", unload=TRUE)
detach("package:ggplot2", unload=TRUE)
detach("package:dplyr", unload=TRUE)
detach("package:cowplot", unload=TRUE)
detach("package:vegan", unload=TRUE)


rm(list=ls(all=TRUE))

sessionInfo()
