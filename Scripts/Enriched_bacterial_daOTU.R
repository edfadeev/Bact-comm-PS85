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

BAC_pruned_DCM <- subset_samples(BAC_pruned, Type =="DCM")

#seperate the fractions
BAC_FL <- subset_samples(BAC_pruned_DCM, Fraction == "0.22")
BAC_FL <- prune_taxa(taxa_sums(BAC_FL)>0,BAC_FL)
BAC_PA <- subset_samples(BAC_pruned_DCM, Fraction == "3")
BAC_PA <- prune_taxa(taxa_sums(BAC_PA)>0,BAC_PA)

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
#Extract daOTU for each region from the DESeq2 analysis (for networks and supp figures)
#####################################
#extract only significant OTU
BAC_FL.DEseq.res.sig <- BAC_FL.DEseq.res[which(BAC_FL.DEseq.res$padj < 0.05), ]
BAC_FL.DEseq.res.sig <- cbind(as(BAC_FL.DEseq.res.sig, "data.frame"),
                              as(tax_table(BAC_FL)[rownames(BAC_FL.DEseq.res.sig), ], "matrix"))
BAC_FL.DEseq.res.sig$Order <- gsub("_unclassified|Incertae Sedis", "_u",BAC_FL.DEseq.res.sig$Order)

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
