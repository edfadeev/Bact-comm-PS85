#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#load libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("plyr"); packageVersion("plyr")
library("dplyr"); packageVersion("dplyr")
library("olsrr"); packageVersion("olsrr")
library("cowplot"); packageVersion("cowplot")
library("iNEXT"); packageVersion("iNEXT")

#set plots theme
theme_set(theme_classic())

#load colour palettes
source('./Scripts/color_palettes.R')

#####################################
#Load phyloseq object
####################################
PS85_BAC <-  readRDS("./Data/BAC_raw.rds")

#####################################
#Alpha diversity
####################################
iNEXT.out <- iNEXT(as.data.frame(otu_table(PS85_BAC)), q=c(0),
                   datatype="abundance", conf = 0.95, nboot = 100)

meta <- as(sample_data(PS85_BAC), "data.frame")
rare <-fortify(iNEXT.out, type=1)
meta$site <- rownames(meta)
rare$Fraction <- meta$Fraction[match(rare$site, meta$site)] 
rare.point <- rare[which(rare$method == "observed"),]
rare.line <- rare[which(rare$method != "observed"),]
rare.line$method <- factor (rare.line$method,
                            c("interpolated", "extrapolated"),
                            c("interpolation", "extrapolation"))

rare.p <- ggplot(rare, aes(x=x, y=y, colour = site))+
  geom_line(aes(linetype = method), lwd = 0.5, data= rare.line)+
  #geom_ribbon(data = rare.line, aes(ymin=y.lwr, ymax= y.upr), alpha = 0.1)+
  geom_point(aes(shape=Fraction), size =3, data= rare.point, colour = "black")+
  scale_colour_discrete(guide = FALSE)+
  scale_x_continuous(limits= c(0,1e+5))+
  labs(x = "Sample size", y = "Species richness")+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

#coverage
df <-fortify(iNEXT.out, type=2)
meta$site <- rownames(meta)
df$Fraction <- meta$Fraction[match(df$site, meta$site)] 

df.point <- df[which(df$method == "observed"),]
df.line <- df[which(df$method != "observed"),]
df.line$method <- factor (df.line$method,
                          c("interpolated", "extrapolated"),
                          c("interpolation", "extrapolation"))


cov.p <- ggplot(df, aes(x=x, y=y, colour = site))+
  geom_line(aes(linetype = method), lwd = 0.5, data= df.line)+
  #geom_ribbon(aes(ymin=y.lwr, ymax= y.upr, colour = NULL), alpha = 0.2)+
  geom_point(aes(shape=Fraction), size =3, data= df.point, colour = "black")+
  scale_colour_discrete(guide = FALSE)+
  scale_x_continuous(limits= c(0,1e+5))+
  scale_y_continuous(breaks=seq(0.9,1,0.05), limits = c(0.9,1))+
  labs(x = "Sample size", y = "Sample coverage")+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

#combined plot
plot_grid(rare.p, cov.p, labels = c("A", "B"), ncol = 2, align = "h")


#####################################
#Table
####################################
BAC_richness <- iNEXT.out$AsyEst[iNEXT.out$AsyEst$Diversity == "Species richness",]
BAC_shannon <- iNEXT.out$AsyEst[iNEXT.out$AsyEst$Diversity == "Shannon diversity",]
BAC_simpson <- iNEXT.out$AsyEst[iNEXT.out$AsyEst$Diversity == "Simpson diversity",]

PS85_comm.char<- data.frame(SampleID = sample_names(PS85_BAC),
                              StationName = sample_data(PS85_BAC)$StationName,
                              Type = sample_data(PS85_BAC)$Type,
                              Fraction = sample_data(PS85_BAC)$Fraction,
                              Sample_sum = iNEXT.out$DataInfo$n,
                              Observed = iNEXT.out$DataInfo$S.obs,
                              Richness = BAC_richness$Observed,
                              Richness.cov = BAC_richness$Observed/BAC_richness$Estimator,
                              Shannon = BAC_shannon$Observed,
                              Shannon.est = BAC_shannon$Estimator,
                              Simpson = BAC_simpson$Observed,
                              Simpson.est = BAC_simpson$Estimator,
                              Sam.comp = 100*iNEXT.out$DataInfo$SC)


write.table(PS85_comm.char, file = "./Data/PS85_comm.char.txt")

#Eukaryotes
PS85_EUK <-  readRDS("./Data/EUK_raw.rds")

#####################################
#Alpha diversity
####################################
iNEXT.out.EUK <- iNEXT(as.data.frame(otu_table(PS85_EUK)), q=c(0), datatype="abundance", conf = 0.95, nboot = 100)

meta <- as(sample_data(PS85_EUK), "data.frame")

rare <-fortify(iNEXT.out.EUK, type=1)
meta$site <- rownames(meta)
rare.point <- rare[which(rare$method == "observed"),]
rare.line <- rare[which(rare$method != "observed"),]
rare.line$method <- factor (rare.line$method,
                            c("interpolated", "extrapolated"),
                            c("interpolation", "extrapolation"))



rare.p <- ggplot(rare, aes(x=x, y=y, colour = site))+
  geom_line(aes(linetype = method), lwd = 0.5, data= rare.line)+
  #geom_ribbon(aes(ymin=y.lwr, ymax= y.upr, colour = NULL), alpha = 0.2)+
  geom_point(size =3, data= rare.point,shape = 3, colour = "black")+
  scale_colour_discrete(guide = FALSE)+
  scale_x_continuous(limits= c(0,1e+5))+
  labs(x = "Sample size", y = "Species richness")+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

#coverage
df <-fortify(iNEXT.out.EUK, type=2)
meta$site <- rownames(meta)

df.point <- df[which(df$method == "observed"),]
df.line <- df[which(df$method != "observed"),]
df.line$method <- factor (df.line$method,
                          c("interpolated", "extrapolated"),
                          c("interpolation", "extrapolation"))


cov.p <- ggplot(df, aes(x=x, y=y, colour = site))+
  geom_line(aes(linetype = method), lwd = 0.5, data= df.line)+
  #geom_ribbon(aes(ymin=y.lwr, ymax= y.upr, colour = NULL), alpha = 0.2)+
  geom_point(size =3, data= df.point, shape = 3, colour = "black")+
  scale_colour_discrete(guide = FALSE)+
  scale_x_continuous(limits= c(0,1e+5))+
  scale_y_continuous(breaks=seq(0.9,1,0.05), limits = c(0.9,1))+
  labs(x = "Sample size", y = "Sample coverage")+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

#combined plot
plot_grid(rare.p, cov.p, labels = c("A", "B"), ncol = 2, align = "h")




EUK_richness <- iNEXT.out.EUK$AsyEst[iNEXT.out.EUK$AsyEst$Diversity == "Species richness",]
EUK_shannon <- iNEXT.out.EUK$AsyEst[iNEXT.out.EUK$AsyEst$Diversity == "Shannon diversity",]



PS85_comm.char.EUK<- data.frame(SampleID = sample_names(PS85_EUK),
                            StationName = sample_data(PS85_EUK)$StationName,
                            Type = sample_data(PS85_EUK)$Type,
                            #Fraction = sample_data(PS85_EUK)$Fraction,
                            Sample_sum = iNEXT.out.EUK$DataInfo$n,
                            Observed = iNEXT.out.EUK$DataInfo$S.obs,
                            Richness = EUK_richness$Estimator,
                            Richness.cov = EUK_richness$Observed/EUK_richness$Estimator,
                            Shannon = EUK_shannon$Estimator,
                            Shannon.cov = EUK_shannon$Observed/EUK_shannon$Estimator,
                            Sam.comp = 100*iNEXT.out.EUK$DataInfo$SC)

write.table(PS85_comm.char.EUK, file = "./Data/PS85_comm-char-EUK.txt")
