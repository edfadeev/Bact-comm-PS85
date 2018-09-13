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
PS85_BAC <-  readRDS("./Data/BAC_pruned.rds")

#####################################
#Alpha diversity
####################################
iNEXT.out <- iNEXT(as.data.frame(otu_table(PS85_BAC)), q=0, datatype="abundance", conf = 0.95, nboot = 100)
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
  #geom_ribbon(aes(ymin=y.lwr, ymax= y.upr, colour = NULL), alpha = 0.2)+
  geom_point(aes(shape=Fraction), size =3, data= rare.point)+
  scale_colour_discrete(guide = FALSE)+
  labs(x = "Sample size", y = "Species richness")+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

#coverage
df <-fortify(iNEXT.out, type=3)
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
  geom_point(aes(shape=Fraction), size =3, data= df.point)+
  scale_colour_discrete(guide = FALSE)+
  scale_x_continuous(breaks=seq(0,1,0.1))+
  labs(x = "Sample coverage", y = "Species richness")+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

#combined plot
plot_grid(rare.p, cov.p, labels = c("A", "B"), ncol = 2, align = "h")

#summary table
PS85_mm <- data.frame()
PS85_alpha_chao1 <- data.frame()
PS85_alpha_chao1.cov <- data.frame()

for (sample in sample_names(PS85_BAC)){
  mm.fit <- rare.line %>%
    filter(method == "interpolation",
           site == sample)
  m1 <- drm(y ~ x, data = mm.fit, fct = MM.2())
  PS85_mm <- rbind(PS85_mm, data.frame(sample, coef(m1)[1]))
 }

PS85_alpha.div <- estimate_richness(PS85_BAC, measures = c("Observed", "Chao1"))

PS85_comm.char<- data.frame(SampleID = sample_names(PS85_BAC),
                              StationName = sample_data(PS85_BAC)$StationName,
                              Type = sample_data(PS85_BAC)$Type,
                              Fraction = sample_data(PS85_BAC)$Fraction,
                              Sample_sum = sample_sums(PS85_BAC),
                              Observed = PS85_alpha.div$Observed,
                              MM.fit = PS85_mm$coef.m1..1.,
                              MM.fit.cov = 100*(PS85_alpha.div$Observed/PS85_mm$coef.m1..1.),
                              Chao = iNEXT.out$AsyEst$Estimator[iNEXT.out$AsyEst$Diversity == "Species richness"],
                              Chao.cov = PS85_alpha.div$Observed/(iNEXT.out$AsyEst$Estimator[iNEXT.out$AsyEst$Diversity == "Species richness"]),
                              Sample.cov = 100*iNEXT.out$DataInfo$SC)


#Eukaryotes
PS85_EUK <-  readRDS("./Data/EUK_pruned.rds")

#####################################
#Alpha diversity
####################################
iNEXT.out <- iNEXT(as.data.frame(otu_table(PS85_EUK)), q=0, datatype="abundance", conf = 0.95, nboot = 100)
meta <- as(sample_data(PS85_EUK), "data.frame")

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
  #geom_ribbon(aes(ymin=y.lwr, ymax= y.upr, colour = NULL), alpha = 0.2)+
  geom_point(aes(shape=Fraction), size =3, data= rare.point)+
  scale_colour_discrete(guide = FALSE)+
  labs(x = "Sample size", y = "Species richness")+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

#coverage
df <-fortify(iNEXT.out, type=3)
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
  geom_point(aes(shape=Fraction), size =3, data= df.point)+
  scale_colour_discrete(guide = FALSE)+
  scale_x_continuous(breaks=seq(0,1,0.1))+
  labs(x = "Sample coverage", y = "Species richness")+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

#combined plot
plot_grid(rare.p, cov.p, labels = c("A", "B"), ncol = 2, align = "h")

#summary table
PS85_mm <- data.frame()
PS85_alpha_chao1 <- data.frame()
PS85_alpha_chao1.cov <- data.frame()

for (sample in sample_names(PS85_EUK)){
  mm.fit <- rare.line %>%
    filter(method == "interpolation",
           site == sample)
  m1 <- drm(y ~ x, data = mm.fit, fct = MM.2())
  PS85_mm <- rbind(PS85_mm, data.frame(sample, coef(m1)[1]))
}

PS85_alpha.div <- estimate_richness(PS85_EUK, measures = c("Observed", "Chao1"))

PS85_comm.char<- data.frame(SampleID = sample_names(PS85_EUK),
                            StationName = sample_data(PS85_EUK)$StationName,
                            Type = sample_data(PS85_EUK)$Type,
                            #Fraction = sample_data(PS85_EUK)$Fraction,
                            Sample_sum = sample_sums(PS85_EUK),
                            Observed = PS85_alpha.div$Observed,
                            MM.fit = PS85_mm$coef.m1..1.,
                            MM.fit.cov = 100*(PS85_alpha.div$Observed/PS85_mm$coef.m1..1.),
                            Chao = iNEXT.out$AsyEst$Estimator[iNEXT.out$AsyEst$Diversity == "Species richness"],
                            Chao.cov = PS85_alpha.div$Observed/(iNEXT.out$AsyEst$Estimator[iNEXT.out$AsyEst$Diversity == "Species richness"]),
                            Sample.cov = 100*iNEXT.out$DataInfo$SC)