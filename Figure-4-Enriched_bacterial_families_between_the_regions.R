#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#load libraries
library("ggplot2"); packageVersion("ggplot2")
library("DESeq2"); packageVersion("DESeq2")
library("phyloseq"); packageVersion("phyloseq")
library("cowplot"); packageVersion("cowplot")
library("gage"); packageVersion("gage")
library("dplyr"); packageVersion("dplyr")

#set plots theme
theme_set(theme_classic())

#load colour palettes
source('./Scripts/color_palettes.R')

#####################################
#Figure 5: Differences in bacterial community composition between the regions
####################################
BAC_pruned<- readRDS("./Data/BAC_pruned.rds")

#seperate the fractions
BAC_FL <- subset_samples(BAC_pruned, Fraction == "0.22")
BAC_FL <- prune_taxa(taxa_sums(BAC_FL)>0,BAC_FL)
BAC_PA <- subset_samples(BAC_pruned, Fraction == "3")
BAC_PA <- prune_taxa(taxa_sums(BAC_PA)>0,BAC_PA)


#define alpha
alpha_val <- 0.05
q_val <- 0.05

##run DEseq2
#run DEseq on FL fraction
BAC_FL.ddsMat <- phyloseq_to_deseq2(BAC_FL, ~Region)
varianceStabilizingTransformation(BAC_FL.ddsMat, blind = TRUE, fitType = "parametric")
BAC_FL.ddsMat <- estimateSizeFactors(BAC_FL.ddsMat)
BAC_FL.ddsMat <- estimateDispersions(BAC_FL.ddsMat)
#geoMeans.FL <- apply(counts(ddsMat_FL), 1, gm_mean)
#ddsMat_FL <- estimateSizeFactors(ddsMat_FL, geoMeans = geoMeans.FL)
BAC_FL.DEseq <- DESeq(BAC_FL.ddsMat, fitType="parametric")
BAC_FL.DEseq.res <- results(BAC_FL.DEseq)

#run DEseq on PA fraction
BAC_PA.ddsMat <- phyloseq_to_deseq2(BAC_PA, ~Region)
varianceStabilizingTransformation(BAC_PA.ddsMat, blind = TRUE, fitType = "parametric")
BAC_PA.ddsMat <- estimateSizeFactors(BAC_PA.ddsMat)
BAC_PA.ddsMat <- estimateDispersions(BAC_PA.ddsMat)
#geoMeans.PA <- apply(counts(ddsMat_PA), 1, gm_mean)
#ddsMat_PA <- estimateSizeFactors(ddsMat_PA, geoMeans = geoMeans.PA)
BAC_PA.DEseq <- DESeq(BAC_PA.ddsMat, fitType="parametric")
BAC_PA.DEseq.res <- results(BAC_PA.DEseq)

#####################################
#create taxonomy db
#####################################
# Identify enriched taxa
#create taxonomy db
BAC_tax <- as.data.frame(tax_table(BAC_pruned))
BAC_tax$OTU <- rownames(tax_table(BAC_pruned))

BAC_tax %>%
  mutate_if(is.factor, as.character) -> BAC_tax

#Extract the desired taxonomic level
Genera <- unique(BAC_tax$Family)
colnames(BAC_tax)[5] <- "id"

#generate list of OTU for each Taxa
OTU.gs <- list()

for (s in 1:length(Genera)){
  n <- Genera[s]
  Order_OTU <- subset(BAC_tax, id == n)
  OTU.gs[[n]] <- Order_OTU$OTU
  
}


#####################################
#Identify enriched families
#####################################

#FL
deseq2.fc <- BAC_FL.DEseq.res$log2FoldChange
names(deseq2.fc) <- rownames(BAC_FL.DEseq.res)
exp.fc=deseq2.fc

fc.Order.p <- gage(exp.fc, gsets = OTU.gs, same.dir=TRUE, ref = NULL, samp = NULL)

#plot
FL_enrch <- rbind(fc.Order.p$greater,fc.Order.p$less)
FL_enrch <-  data.frame(FL_enrch[FL_enrch[,"q.val"]<0.05 &
                                   !is.na(FL_enrch[,"q.val"]),])
FL_enrch$id <- rownames(FL_enrch)
FL_enrch <- unique(merge(FL_enrch,BAC_tax[,c("id","Class")]))

#PA
deseq2.fc_PA <- BAC_PA.DEseq.res$log2FoldChange
names(deseq2.fc_PA) <- rownames(BAC_PA.DEseq.res)
exp.fc_PA=deseq2.fc_PA


fc.Order.p_PA <- gage(exp.fc_PA, gsets = OTU.gs, same.dir=TRUE, ref = NULL, samp = NULL,compare = "unpaired")

#plot
PA_enrch <- rbind(fc.Order.p_PA$greater,fc.Order.p_PA$less)
PA_enrch <-  data.frame(PA_enrch[PA_enrch[,"q.val"]<0.05 &
                                   !is.na(PA_enrch[,"q.val"]),])
PA_enrch$id <- rownames(PA_enrch)
PA_enrch <- unique(merge(PA_enrch,BAC_tax[,c("id","Class")]))


#merge dataframes and plot
FL_enrch$frac <- "FL"
PA_enrch$frac <- "PA"

both_enr <- rbind(FL_enrch,PA_enrch)

#correct taxa
both_enr$id <- gsub("_unclassified", "_unc",both_enr$id)
both_enr$id <- gsub("(Marine group B)", "",both_enr$id)
both_enr$id <- gsub("(SAR406\ clade)", "",both_enr$id)
both_enr$Class <- gsub("_unclassified", "_unc",both_enr$Class)

#plot
enrch.plot <- ggplot(data=both_enr, aes(y=stat.mean , x=id, label = set.size))+ 
  geom_text(size = 4, aes(y=stat.mean , x=id), nudge_y= -1.4, nudge_x= 0)+
  ylab("Mean log2foldchange")+ 
  geom_point(size = 5, aes(colour = Class))+
  scale_colour_manual(values = tol21rainbow)+ 
  ylim(-10,10)+
  scale_x_discrete("Family",expand = waiver())+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  theme(legend.position = "bottom")+
  facet_wrap(~frac)+
  coord_flip()

ggsave("./Figures/Figure-4-enriched_taxa_bac.pdf", enrch.plot,
       dpi = 300, device= "pdf",width = 17, height = 22, units = "cm")

#####################################
#explore enriched taxa
#####################################
# Get the 5 most enriched in WSC-FL
BAC_FL.enr.WSC <- data.frame(id=rownames(fc.Order.p$greater), fc.Order.p$greater) %>% 
  tbl_df() %>% 
  dplyr::filter(q.val<0.05,row_number()<=5) %>% 
  .$id %>% 
  as.character()

# Get the 5 most enriched in EGC-FL
BAC_FL.enr.EGC <- data.frame(id=rownames(fc.Order.p$less), fc.Order.p$less) %>% 
  tbl_df() %>% 
  dplyr::filter(q.val<0.05,row_number()<=5) %>% 
  .$id %>% 
  as.character()

# Get the 5 most enriched in WSC-PA
BAC_PA.enr.WSC <- data.frame(id=rownames(fc.Order.p_PA$greater), fc.Order.p_PA$greater) %>% 
  tbl_df() %>% 
  dplyr::filter(q.val<0.05,row_number()<=5) %>% 
  .$id %>% 
  as.character()

# Get the 5 most enriched in EGC-PA
BAC_PA.enr.EGC <- data.frame(id=rownames(fc.Order.p_PA$less), fc.Order.p_PA$less) %>% 
  tbl_df() %>% 
  dplyr::filter(q.val<0.05,row_number()<=5) %>% 
  .$id %>% 
  as.character()



#####################################
#Extract daOTU for each region from the DESeq2 analysis (for networks and supp figures)
#####################################
#extract only significant OTU
BAC_FL.DEseq.res.sig <- BAC_FL.DEseq.res[which(BAC_FL.DEseq.res$padj < 0.05), ]
BAC_FL.DEseq.res.sig <- cbind(as(BAC_FL.DEseq.res.sig, "data.frame"),
                              as(tax_table(BAC_FL)[rownames(BAC_FL.DEseq.res.sig), ], "matrix"))
BAC_FL.DEseq.res.sig$Order <- gsub("_unclassified| Incertae Sedis", "_unc",BAC_FL.DEseq.res.sig$Order)

#devide the OTUs by the region they were enriched in
BAC_FL.DEseq.WSC  <- BAC_FL.DEseq.res.sig[BAC_FL.DEseq.res.sig[, "log2FoldChange"] > 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Genus") ]
BAC_FL.DEseq.EGC  <- BAC_FL.DEseq.res.sig[BAC_FL.DEseq.res.sig[, "log2FoldChange"] < 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Genus")]

#extract only significant OTU
BAC_PA.DEseq.res.sig <- BAC_PA.DEseq.res[which(BAC_PA.DEseq.res$padj < 0.05), ]
BAC_PA.DEseq.res.sig <- cbind(as(BAC_PA.DEseq.res.sig, "data.frame"),
                              as(tax_table(BAC_PA)[rownames(BAC_PA.DEseq.res.sig), ], "matrix"))
BAC_PA.DEseq.res.sig$Order <- gsub("_unclassified| Incertae Sedis", "_unc",BAC_PA.DEseq.res.sig$Order)

#devide the OTUs by the region they were enriched in
BAC_PA.DEseq.WSC  <- BAC_PA.DEseq.res.sig[BAC_PA.DEseq.res.sig[, "log2FoldChange"] > 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Genus") ]
BAC_PA.DEseq.EGC  <- BAC_PA.DEseq.res.sig[BAC_PA.DEseq.res.sig[, "log2FoldChange"] < 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Genus")]

#list the enriched OTUs from DESeq analysis
FL.WSC.enr <- data.frame(name=rownames(BAC_FL.DEseq.WSC), region = "FL.WSC")
PA.WSC.enr <- data.frame(name=rownames(BAC_PA.DEseq.WSC), region = "PA.WSC")

FL.EGC.enr <- data.frame(name=rownames(BAC_FL.DEseq.EGC), region = "FL.EGC")
PA.EGC.enr <- data.frame(name=rownames(BAC_PA.DEseq.EGC), region = "PA.EGC")                  

BAC_daOTU <- unique(rbind(FL.WSC.enr,PA.WSC.enr,FL.EGC.enr,PA.EGC.enr))


#####################################
#Supplementary Figure 5. Differences in bacterial community composition between the regions
#####################################
#Aggregate on Order level and bind both regions
BAC_FL.DEseq.WSC.agg <-as.data.frame(as.list(aggregate(log2FoldChange~Class+Order, 
                                                       BAC_FL.DEseq.WSC, 
                                                       FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))
BAC_FL.DEseq.EGC.agg <-as.data.frame(as.list(aggregate(log2FoldChange~Class+Order, 
                                                       BAC_FL.DEseq.EGC, 
                                                       FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))

FL_enriched_agg <- rbind(BAC_FL.DEseq.WSC.agg,BAC_FL.DEseq.EGC.agg)
FL_enriched_agg %>% mutate_if(is.factor, as.character) -> FL_enriched_agg

BAC_PA.DEseq.WSC.agg <-as.data.frame(as.list(aggregate(log2FoldChange~Class+Order, 
                                                       BAC_PA.DEseq.WSC, 
                                                       FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))
BAC_PA.DEseq.EGC.agg <-as.data.frame(as.list(aggregate(log2FoldChange~Class+Order, 
                                                       BAC_PA.DEseq.EGC, 
                                                       FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))
PA_enriched_agg <- rbind(BAC_PA.DEseq.WSC.agg,BAC_PA.DEseq.EGC.agg)
PA_enriched_agg %>% mutate_if(is.factor, as.character) -> PA_enriched_agg

#plot only Orders with at least 5 enriched OTU
BAC_FL.DEseq.agg <-as.data.frame(as.list(aggregate(log2FoldChange~Class+Order, 
                                                   BAC_FL.DEseq.res.sig, FUN = function(x) count=length(x))))
BAC_FL.DEseq.agg <- BAC_FL.DEseq.agg[BAC_FL.DEseq.agg$log2FoldChange>10,]
FL_enriched_agg <- FL_enriched_agg[FL_enriched_agg$Order %in% BAC_FL.DEseq.agg$Order,]


BAC_PA.DEseq.agg <-as.data.frame(as.list(aggregate(log2FoldChange~Class+Order, 
                                                   BAC_PA.DEseq.res.sig, FUN = function(x) count=length(x))))
BAC_PA.DEseq.agg <- BAC_PA.DEseq.agg[BAC_PA.DEseq.agg$log2FoldChange>10,]
PA_enriched_agg <- PA_enriched_agg[PA_enriched_agg$Order %in% BAC_PA.DEseq.agg$Order,]

#merged the fractions
FL_enriched_agg$frac <- "FL"
PA_enriched_agg$frac <- "PA"

BAC_daOTU.agg <- rbind(FL_enriched_agg,PA_enriched_agg)

#plot
BAC_daOTU.p <- ggplot(data=BAC_daOTU.agg, aes(y=log2FoldChange.mean , x=Order, fill = Class, label = log2FoldChange.count))+ 
  geom_text(aes(y=log2FoldChange.mean , x=Order), nudge_y= 0, nudge_x= -0.45)+
  geom_errorbar(aes(ymin = log2FoldChange.mean-log2FoldChange.se, ymax = log2FoldChange.mean +log2FoldChange.se), width = 0.5) +   
  ylab("log2foldchange")+
  geom_point(size = 5, shape = 21)+
  scale_x_discrete("Order")+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  theme(legend.position = "bottom")+
  scale_fill_manual(values = phyla.col)+
  facet_wrap(~frac)+
  coord_flip() 

ggsave("./Figures/Supp-Figure-4-daOTU.png", BAC_daOTU.p,
       dpi = 300, device= "png",width = 30, height = 30, units = "cm")


###explore the results
#taxa overlaps
#overview the results
summary(BAC_FL.DEseq.res)
summary(BAC_PA.DEseq.res)

#OTU overlaps
#extract only OTU with pvalue
BAC_FL.DEseq.res.padj <- BAC_FL.DEseq.res[!is.na(BAC_FL.DEseq.res$padj),]
BAC_PA.DEseq.res.padj <- BAC_PA.DEseq.res[!is.na(BAC_PA.DEseq.res$padj),]


#overlaping
FL_over <- intersect(BAC_FL.DEseq.WSC.agg$Order,BAC_FL.DEseq.EGC.agg$Order)
PA_over <- intersect(BAC_PA.DEseq.WSC.agg$Order,BAC_PA.DEseq.EGC.agg$Order)
#unique
FL_WSC_unique <- setdiff(BAC_FL.DEseq.WSC.agg$Order,BAC_FL.DEseq.EGC.agg$Order)
FL_EGC_unique <- setdiff(BAC_FL.DEseq.EGC.agg$Order,BAC_FL.DEseq.WSC.agg$Order)

PA_WSC_unique <- setdiff(BAC_PA.DEseq.WSC.agg$Order,BAC_PA.DEseq.EGC.agg$Order)
PA_EGC_unique <- setdiff(BAC_PA.DEseq.EGC.agg$Order,BAC_PA.DEseq.WSC.agg$Order)

#identify the order with most enriched OTUs
BAC_FL.DEseq.res.sig.top <- cbind(BAC_FL.DEseq.res.sig, otu = rownames(BAC_FL.DEseq.res.sig))
BAC_PA.DEseq.res.sig.top <- cbind(BAC_PA.DEseq.res.sig, otu = rownames(BAC_PA.DEseq.res.sig))


top.orders.FL <- 
  as.data.frame(as.list(aggregate(otu~Class+Order, 
                                  BAC_FL.DEseq.res.sig.top, 
                                  FUN = function(x) count=length(x))))


top.orders.PA <- 
  as.data.frame(as.list(aggregate(otu~Class+Order, 
                                  BAC_PA.DEseq.res.sig.top, 
                                  FUN = function(x) count=length(x))))

#####################################
#Supplementary Figure 5: Enriched microbial eukaryotic taxonomic groups between the regions
#####################################
EUK_pruned<- readRDS("./Data/EUK_pruned.rds")

#run DEseq on EUK fraction
EUK_ddsMat <- phyloseq_to_deseq2(EUK_pruned, ~Region)
varianceStabilizingTransformation(EUK_ddsMat, blind = TRUE, fitType = "parametric")
EUK_ddsMat <- estimateSizeFactors(EUK_ddsMat)
EUK_ddsMat <- estimateDispersions(EUK_ddsMat)
EUK_DEseq <- DESeq(EUK_ddsMat, fitType="parametric")
EUK_DEseq.res <- results(EUK_DEseq)

#####################################
# Identify enriched taxa
#####################################
# Identify enriched taxa
#create taxonomy db
EUK_tax <- as.data.frame(tax_table(EUK_pruned))
EUK_tax$OTU <- rownames(tax_table(EUK_pruned))

EUK_tax %>%
  mutate_if(is.factor, as.character) -> EUK_tax

#Extract the desired taxonomic level
Genera <- unique(EUK_tax$X5)
colnames(EUK_tax)[5] <- "id"

#generate list of OTU for each Taxa
OTU.gs_EUK <- list()

for (s in 1:length(Genera)){
  n <- Genera[s]
  Order_OTU <- subset(EUK_tax, id == n)
  OTU.gs_EUK[[n]] <- Order_OTU$OTU
  
}

#####################################
#Identify enriched families
#####################################
EUK_deseq2.fc <- EUK_DEseq.res$log2FoldChange
names(EUK_deseq2.fc) <- rownames(EUK_DEseq.res)

EUK_fc.Order.p <- gage(EUK_deseq2.fc, gsets = OTU.gs_EUK, same.dir=TRUE, ref = NULL, samp = NULL)

#plot
EUK_enrch <- rbind(EUK_fc.Order.p$greater,EUK_fc.Order.p$less)
EUK_enrch.sig <-  data.frame(EUK_enrch[EUK_enrch[,"q.val"]<0.05 &
                                   !is.na(EUK_enrch[,"q.val"]),])
EUK_enrch.sig$id <- rownames(EUK_enrch.sig)

EUK_enr.taxa.sig<- unique(merge(EUK_enrch.sig, EUK_tax[,c("id","X4")]))

EUK_enr.taxa.sig.plot <- ggplot(data=EUK_enr.taxa.sig, aes(y=stat.mean , x=id, fill = X4, label = set.size))+ 
  geom_text(aes(y=stat.mean , x=id), nudge_y= 0, nudge_x= -0.45)+
  ylab("log2foldchange")+
  geom_point(size = 5, shape = 21)+
  #ylim(-12,12)+
  scale_x_discrete("Taxa")+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  theme(legend.position = "bottom")+
  #scale_fill_manual(values = phyla.col)+
  coord_flip() 

#####################################
#Extract daOTU for each region from the DESeq2 analysis (for networks and supp figures)
#####################################
#extract only significant OTU
EUK.DEseq.res.sig <- EUK_DEseq.res[which(EUK_DEseq.res$padj < 0.05), ]
EUK.DEseq.res.sig <- cbind(as(EUK.DEseq.res.sig, "data.frame"),
                              as(tax_table(EUK_pruned)[rownames(EUK.DEseq.res.sig), ], "matrix"))
EUK.DEseq.res.sig$X3 <- gsub("_unclassified| Incertae Sedis", "_unc",EUK.DEseq.res.sig$X3)

#devide the OTUs by the region they were enriched in
EUK.DEseq.WSC  <- EUK.DEseq.res.sig[EUK.DEseq.res.sig[, "log2FoldChange"] > 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "X2", "X3", "X4", "X5", "X6") ]
EUK.DEseq.EGC  <- EUK.DEseq.res.sig[EUK.DEseq.res.sig[, "log2FoldChange"] < 0,c("baseMean", "log2FoldChange", "lfcSE", "padj",  "X2", "X3", "X4", "X5", "X6")]


#list the enriched OTUs from DESeq analysis
EUK.WSC.enr <- data.frame(name=rownames(EUK.DEseq.WSC), region = "WSC")
EUK.EGC.enr <- data.frame(name=rownames(EUK.DEseq.EGC), region = "EGC")                  

EUK_daOTU <- unique(rbind(EUK.WSC.enr,EUK.EGC.enr))


#Aggregate on Order level and bind both regions
EUK.DEseq.WSC.agg <-as.data.frame(as.list(aggregate(log2FoldChange~X3+X4+X5, 
                                                       EUK.DEseq.WSC, 
                                                       FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))
EUK.DEseq.EGC.agg <-as.data.frame(as.list(aggregate(log2FoldChange~X3+X4+X5, 
                                                       EUK.DEseq.EGC, 
                                                       FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))

FL_enriched_agg <- rbind(EUK.DEseq.WSC.agg,EUK.DEseq.EGC.agg)
FL_enriched_agg %>% mutate_if(is.factor, as.character) -> FL_enriched_agg


#plot only Orders with at least 5 enriched OTU
EUK.DEseq.agg <-as.data.frame(as.list(aggregate(log2FoldChange~X3+X4+X5, 
                                                   EUK.DEseq.res.sig, FUN = function(x) count=length(x))))
EUK.DEseq.agg <- EUK.DEseq.agg[EUK.DEseq.agg$log2FoldChange>3,]
FL_enriched_agg <- FL_enriched_agg[FL_enriched_agg$X5 %in% EUK.DEseq.agg$X5,]

#plot
EUK_daOTU.p <- ggplot(data=FL_enriched_agg, aes(y=log2FoldChange.mean , x=X5, fill = X3, label = log2FoldChange.count))+ 
  geom_text(aes(y=log2FoldChange.mean , x=X5), nudge_y= 0, nudge_x= -0.45)+
  geom_errorbar(aes(ymin = log2FoldChange.mean-log2FoldChange.se, ymax = log2FoldChange.mean +log2FoldChange.se), width = 0.5) +   
  ylab("log2foldchange")+
  geom_point(size = 5, shape = 21)+
  scale_x_discrete("X5")+ ylim(-9,9)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  theme(legend.position = "bottom")+
  #scale_fill_manual(values = phyla.col)+
  coord_flip() 

ggsave("./Figures/Supp-Figure-5-Euk_daOTU.png", EUK_daOTU.p,
       dpi = 300, device= "png",width = 30, height = 30, units = "cm")

#####################################
#remove all temporary datasets 
####################################
detach("package:phyloseq", unload=TRUE)
detach("package:UpSetR", unload=TRUE)
detach("package:ggplot2", unload=TRUE)
detach("package:DESeq2", unload=TRUE)
detach("package:cowplot", unload=TRUE)
detach("package:gage", unload=TRUE)
detach("package:dplyr", unload=TRUE)

rm(list=ls(all=TRUE))

sessionInfo()
