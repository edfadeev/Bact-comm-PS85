#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#install required packages
.inst_packages <- c("ggplot2", "gridExtra", "phyloseq",
                    "igraph","UpSetR","vegan","cowplot",
                    "dplyr","ggpmisc","VennDiagram",
                    "gage","DESeq2","zoo")

.inst <- .inst_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.inst_packages, ask = F)
}

#load libraries
library("phyloseq"); packageVersion("phyloseq")
library("DESeq2"); packageVersion("DESeq2")
library("zoo"); packageVersion("zoo")

#load functions
source('./Scripts/ReadAmplicon.R')
source('./Scripts/ReadAmplicon_EUK.R')

#####################################
#Parse bacterial dataset for Phyloseq
#####################################
BAC <- ReadAmplicon(otu = "./Data/Bacteria/OTU_contingency_table.csv", tax = "./Data/Bacteria/amplicons_seeds_taxonomy.txt", 
                    silva = "./Data/Silva_tax/SILVA128_tax_slv_ssu_curated.tsv", 
                    domain = "Bacteria", silva.sep = "\t", singletons = F, unclassified = T, write.files = F)
OTU_BAC<- BAC$OTU
TAX_BAC<- BAC$TAX

# Contextual data
ENV_BAC <- read.csv("./Data/Bacteria/sample_list_PS85_bact_subset.csv", sep = "," , h = T, row.names = 1, fill = T, na.strings=c("","NA"))

# Check order of samples
all.equal(rownames(OTU_BAC), rownames(TAX_BAC))

#creating Phyloseq dataset
BAC_raw <- phyloseq(otu_table(OTU_BAC, taxa_are_rows = TRUE), tax_table(TAX_BAC), sample_data(ENV_BAC))

#Fix categories 
sample_data(BAC_raw)$Type <- factor(
  sample_data(BAC_raw)$Type, 
  levels = c("SRF", "DCM", "BDCM", "MESO", "UNDEF"))

sample_data(BAC_raw)$Water.mass<- factor(
  sample_data(BAC_raw)$Water.mass, 
  levels = c("AW", "PSW", "PSWw", "#REF!"))

sample_data(BAC_raw)$Fraction<- factor(
  sample_data(BAC_raw)$Fraction, 
  levels = c("0.22", "3"))

sample_data(BAC_raw)$StationName<- factor(
  sample_data(BAC_raw)$StationName, 
  levels = c("10W","8.5W","7W",
             "EG1","EG2","EG3","EG4","1W", "1E", "HG9", 
             "HG4", "HG1","SV4","SV3","SV2","SV1","N4","S3"))

#Fix taxonomy 
colnames(BAC_raw@tax_table)<- c("Domain","Phylum","Class","Order","Family","Genus")
BAC_raw@tax_table[is.na(BAC_raw@tax_table)]<- "Unclassified"

#####################################
#Preprocess bacterial OTU by prevalence of each taxa
#####################################
#in how many samples did each taxa appear at least once
prev0 = apply(X = otu_table(BAC_raw),
              MARGIN = ifelse(taxa_are_rows(BAC_raw), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(BAC_raw)

# Execute prevalence filter, using `prune_taxa()` function
BAC_pruned <-  prune_taxa((prev0 > prevalenceThreshold), BAC_raw)

#######################################
# Variance stabilized Transformation
#######################################
#after prevalance
BAC_pruned.dds <- phyloseq_to_deseq2(BAC_pruned, ~1)
varianceStabilizingTransformation(BAC_pruned.dds, blind = TRUE, fitType = "parametric")
BAC_pruned.dds <- estimateSizeFactors(BAC_pruned.dds)
BAC_pruned.dds <- estimateDispersions(BAC_pruned.dds)
otu.vst <- getVarianceStabilizedData(BAC_pruned.dds)

#make sure that the dimentions aof the OTU table and the DEseq object are matching
dim(otu.vst)
dim(otu_table(BAC_pruned))

BAC_pruned.vst<-BAC_pruned
otu_table(BAC_pruned.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)

saveRDS(BAC_pruned.vst, "./Data/BAC_pruned_vst.rds")
saveRDS(BAC_pruned, "./Data/BAC_pruned.rds")
saveRDS(BAC_raw, "./Data/BAC_raw.rds")

#####################################
#Parse eukaryotic dataset for Phyloseq
#####################################
EUK <- ReadAmplicon_EUK(otu = "./Data/Eukaryotes/OTU_contingency_table.csv", 
                        tax = "./Data/Eukaryotes/amplicons_seeds_taxonomy_curated.txt", 
                        silva = "./Data/Silva_tax/SILVA128_tax_slv_EUK_curated.txt", 
                        domain = "Eukaryota", silva.sep = "\t", singletons = F, unclassified = T, write.files = F)

#adjust manual taxonomy
TAX_EUK<- as.data.frame(read.csv("./Data/Eukaryotes/curated_taxonomy_euk.csv",sep = "," , h = T, row.names = 1, strip.white=T))
TAX_EUK <- subset(TAX_EUK, X1!="#N/A")
TAX_EUK <- subset(TAX_EUK, X1!="Unclassified")

TAX_EUK_t <- as.data.frame(t(TAX_EUK))
TAX_EUK_t[TAX_EUK_t=="0"] <- ""
TAX_EUK_t %>% do(zoo::na.locf(.)) -> TAX_EUK_t
TAX_EUK <- t(TAX_EUK_t)

#remove unclassifed OTU
OTU_EUK<- EUK$OTU[c(rownames(TAX_EUK)),]

#Contextual data
ENV_EUK <- read.csv("./Data/Eukaryotes/sample_list_PS85_euk_subset.csv", sep = "," , h = T, row.names = 1, fill = T, na.strings=c("","NA"))


#Check order of samples
all.equal(rownames(OTU_EUK), rownames(TAX_EUK))

#creating Phyloseq dataset
EUK_raw <- phyloseq(otu_table(OTU_EUK, taxa_are_rows = TRUE), 
                    tax_table(as.matrix(TAX_EUK)), 
                    sample_data(ENV_EUK))

#Fix categories
sample_data(EUK_raw)$Type <- factor(
  sample_data(EUK_raw)$Type, 
  levels = c("SRF", "DCM", "BDCM", "MESO", "UNDEF"))

sample_data(EUK_raw)$Water.mass<- factor(
  sample_data(EUK_raw)$Water.mass, 
  levels = c("AW", "PSW", "PSWw", "#REF!"))

sample_data(EUK_raw)$StationName<- factor(
  sample_data(EUK_raw)$StationName, 
  levels = c("10W","8.5W","7W",
             "EG1","EG2","EG3","EG4","1W", "1E", "HG9", 
             "HG4", "HG1","SV4","SV3","SV2","SV1","N4","S3"))

#remove Eumetazoa
EUK_raw <- subset_taxa(EUK_raw, X3!= "Holozoa")

#####################################
#Preprocess eukaryotic OTU by prevalence of each taxa
#####################################
#in how many samples did each taxa appear at least once
prev0_EUK = apply(X = otu_table(EUK_raw),
                  MARGIN = ifelse(taxa_are_rows(EUK_raw), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})

# Define prevalence threshold as 5% of total samples
prevalenceThreshold_EUK = 0.05 * nsamples(EUK_raw)
prevalenceThreshold_EUK
# Execute prevalence filter, using `prune_taxa()` function
EUK_pruned <-  prune_taxa((prev0_EUK > prevalenceThreshold_EUK), EUK_raw)

#######################################
# Variance stabilized Transformation
#######################################
#after prevalance
EUK_pruned.dds <- phyloseq_to_deseq2(EUK_pruned, ~1)
varianceStabilizingTransformation(EUK_pruned.dds, blind = TRUE, fitType = "parametric")
EUK_pruned.dds <- estimateSizeFactors(EUK_pruned.dds)
EUK_pruned.dds <- estimateDispersions(EUK_pruned.dds)
otu.vst <- getVarianceStabilizedData(EUK_pruned.dds)

#make sure that the dimentions aof the OTU table and the DEseq object are matching
dim(otu.vst)
dim(otu_table(EUK_pruned))

EUK_pruned.vst<-EUK_pruned
otu_table(EUK_pruned.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)

saveRDS(EUK_pruned.vst, "./Data/EUK_pruned_vst.rds")
saveRDS(EUK_pruned, "./Data/EUK_pruned.rds")
saveRDS(EUK_raw, "./Data/EUK_raw.rds")
#####################################
#export for network analysis
####################################
#run the 'source' to export  OTU tables and taxonomy for CoNet
source("./Scripts/export_for_conet.R")
#the files are manualy edited in excel (need to add a cell above the row names)
#the OTU tables are feeded into CoNet package in cytoscape
#################################################################

#####################################
#remove all temporary datasets 
####################################
detach("package:phyloseq", unload=TRUE)
rm(list=ls(all=TRUE))

sessionInfo()
