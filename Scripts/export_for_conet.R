#####################################
#export OTU tables for CoNet
#####################################
#extract only those who have bacterial samples with EUK samples
BAC_pruned_EUK<- BAC_pruned %>% # SELECT THE PROPER DATASET
  subset_samples(
    !is.na(EUK_sampaleID))

EUK_phytoplankton_glom_pruned<- prune_samples(levels(sample_data(BAC_pruned_EUK)$EUK_sampaleID),EUK_pruned)

######################################
#FL
######################################
BAC_WSC_FL <- subset_samples(BAC_pruned_EUK, Fraction == "0.22")
BAC_WSC_FL <- prune_taxa(taxa_sums(BAC_WSC_FL)>0, BAC_WSC_FL)
write.table(t(veganotu(BAC_WSC_FL)),file="./Data/network/for_CoNet/BAC_FL_all.txt",sep="\t", quote=FALSE)

EUK_phyto_WSC_FL <- prune_samples(levels(sample_data(BAC_WSC_FL)$EUK_sampaleID), EUK_phytoplankton_glom_pruned)

EUK_phyto_WSC_FL_otu <- as.data.frame(veganotu(EUK_phyto_WSC_FL))

EUK_phyto_WSC_FL_otu$Euk_sampleID <- rownames(EUK_phyto_WSC_FL_otu)

rownames(EUK_phyto_WSC_FL_otu)<- join(EUK_phyto_WSC_FL_otu, data.frame(Euk_sampleID = sample_data(BAC_WSC_FL)$EUK_sampaleID, BAC_sampleID= rownames(sample_data(BAC_WSC_FL))), by = "Euk_sampleID" )[,c("BAC_sampleID")]

EUK_phyto_WSC_FL_otu <- subset(EUK_phyto_WSC_FL_otu, select = -c(Euk_sampleID))

write.table(t(EUK_phyto_WSC_FL_otu),file="./Data/network/for_CoNet/EUK_FL_all.txt",sep="\t", quote=FALSE)

######################################
#PA
######################################
BAC_WSC_PA <- subset_samples(BAC_pruned_EUK, Fraction == "3")
BAC_WSC_PA <- prune_taxa(taxa_sums(BAC_WSC_PA)>0, BAC_WSC_PA)
write.table(t(veganotu(BAC_WSC_PA)),file="./Data/network/for_CoNet/BAC_PA_all.txt",sep="\t", quote=FALSE)


EUK_phyto_WSC_PA <- prune_samples(levels(sample_data(BAC_WSC_PA)$EUK_sampaleID), EUK_phytoplankton_glom_pruned)

EUK_phyto_WSC_PA_otu <- as.data.frame(veganotu(EUK_phyto_WSC_PA))

EUK_phyto_WSC_PA_otu$Euk_sampleID <- rownames(EUK_phyto_WSC_PA_otu)

rownames(EUK_phyto_WSC_PA_otu)<- join(EUK_phyto_WSC_PA_otu, data.frame(Euk_sampleID = sample_data(BAC_WSC_PA)$EUK_sampaleID, BAC_sampleID= rownames(sample_data(BAC_WSC_PA))), by = "Euk_sampleID" )[,c("BAC_sampleID")]

EUK_phyto_WSC_PA_otu <- subset(EUK_phyto_WSC_PA_otu, select = -c(Euk_sampleID))

write.table(t(EUK_phyto_WSC_PA_otu),file="./Data/network/for_CoNet/EUK_PA_all.txt",sep="\t", quote=FALSE)

##############################################################################
#the files are manualy edited in excel 
#the otu tables are feeded into CoNet package in cytoscape
##############################################################################