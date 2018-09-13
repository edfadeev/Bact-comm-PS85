#load libraries
library(dplyr)
library(plyr)

BAC_DCM <- subset_samples(BAC_pruned, Type == "DCM")
EUK_DCM <- subset_samples(EUK_pruned, Type == "DCM")
#####################################
#export OTU tables for CoNet
#####################################
#extract only those who have bacterial samples with EUK samples
BAC_DCM_EUK<- BAC_DCM %>% # SELECT THE PROPER DATASET
  subset_samples(
    !is.na(EUK_sampaleID))

EUK_phytoplankton_glom_DCM<- prune_samples(levels(sample_data(BAC_DCM_EUK)$EUK_sampaleID),EUK_DCM)

######################################
#FL
######################################
BAC_WSC_FL <- subset_samples(BAC_DCM_EUK, Fraction == "0.22")
BAC_WSC_FL <- prune_taxa(taxa_sums(BAC_WSC_FL)>0, BAC_WSC_FL)
write.table(t(veganotu(BAC_WSC_FL)),file="./Data/network/for_CoNet/BAC_FL_DCM.txt",sep="\t", quote=FALSE)

EUK_phyto_WSC_FL <- prune_samples(levels(sample_data(BAC_WSC_FL)$EUK_sampaleID), EUK_phytoplankton_glom_DCM)

EUK_phyto_WSC_FL_otu <- as.data.frame(veganotu(EUK_phyto_WSC_FL))

EUK_phyto_WSC_FL_otu$Euk_sampleID <- rownames(EUK_phyto_WSC_FL_otu)

rownames(EUK_phyto_WSC_FL_otu)<- join(EUK_phyto_WSC_FL_otu, data.frame(Euk_sampleID = sample_data(BAC_WSC_FL)$EUK_sampaleID, BAC_sampleID= rownames(sample_data(BAC_WSC_FL))), by = "Euk_sampleID" )[,c("BAC_sampleID")]

EUK_phyto_WSC_FL_otu <- subset(EUK_phyto_WSC_FL_otu, select = -c(Euk_sampleID))

write.table(t(EUK_phyto_WSC_FL_otu),file="./Data/network/for_CoNet/EUK_FL_DCM.txt",sep="\t", quote=FALSE)

######################################
#PA
######################################
BAC_WSC_PA <- subset_samples(BAC_DCM_EUK, Fraction == "3")
BAC_WSC_PA <- prune_taxa(taxa_sums(BAC_WSC_PA)>0, BAC_WSC_PA)
write.table(t(veganotu(BAC_WSC_PA)),file="./Data/network/for_CoNet/BAC_PA_DCM.txt",sep="\t", quote=FALSE)


EUK_phyto_WSC_PA <- prune_samples(levels(sample_data(BAC_WSC_PA)$EUK_sampaleID), EUK_phytoplankton_glom_DCM)

EUK_phyto_WSC_PA_otu <- as.data.frame(veganotu(EUK_phyto_WSC_PA))

EUK_phyto_WSC_PA_otu$Euk_sampleID <- rownames(EUK_phyto_WSC_PA_otu)

rownames(EUK_phyto_WSC_PA_otu)<- join(EUK_phyto_WSC_PA_otu, data.frame(Euk_sampleID = sample_data(BAC_WSC_PA)$EUK_sampaleID, BAC_sampleID= rownames(sample_data(BAC_WSC_PA))), by = "Euk_sampleID" )[,c("BAC_sampleID")]

EUK_phyto_WSC_PA_otu <- subset(EUK_phyto_WSC_PA_otu, select = -c(Euk_sampleID))

write.table(t(EUK_phyto_WSC_PA_otu),file="./Data/network/for_CoNet/EUK_PA_DCM.txt",sep="\t", quote=FALSE)


########################################Generate taxonomy table#############################
PS85_TAX <- as.data.frame(rbind(tax_table(BAC_DCM), tax_table(EUK_DCM)))
PS85_TAX$OTUs <- rownames(PS85_TAX)
PS85_TAX$Taxon <- rownames(PS85_TAX)
PS85_TAX$Lineage <- paste(PS85_TAX$Domain,PS85_TAX$Phylum,PS85_TAX$Class,
                          PS85_TAX$Order,PS85_TAX$Family,PS85_TAX$Genus,PS85_TAX$OTUs, sep= ";")

PS85_TAX <- PS85_TAX[,c("OTUs","Lineage","Domain","Phylum","Class","Order","Family","Genus")]

write.table(PS85_TAX,file="./Data/network/for_CoNet/PS85_tax_net.txt",sep="\t", quote=FALSE, row.names = F)

##############################################################################
#the files are manualy edited in excel 
#the otu tables are feeded into CoNet package in cytoscape
##############################################################################