This GitHub Repository contains code and data included in:

Fadeev, E., Salter, I., Schourup-Kristensen, V., Nöthig, E.-M., Metfies, K., Engel, A., et al. (2018). **Microbial communities in the East and West Fram Strait during sea-ice melting season**. Front. Mar. Sci. 5, 429. doi:10.3389/fmars.2018.00429.

## Data:
The data directory includes the required otu table for both bacterial and eukaryotic communities, as well as the output from CoNet for generating the subnetworks.

## Files:
* *dataset_preprocess.R* - run at the beginning of the analysis
* *Alpha_diversity_calculations.R* - produces the rarefactions and the alpha diversity indeces of the microbial communities (Table S1 in the paper)
* *Figure-3-OTU_overlaps_between_regions.R* - calculates shared OTUs between microbial communities (Figure 3 in the paper)
* *Figure-5-Enriched_bacterial_families_between_the_regions.R* - calculates using DESeq2 and GAGE differentially abundant OTUs and enriched taxa (Figure 5 in the paper)
* *Figure-7-Bacterial_community_characteristics.R* - calculates dissimilarity between microbial communities (Figure 7 in the paper)
* *Figure-8-RDA_ordination_of_bacterial_community_composition.R* - RDA analysis for identification of the main environmental drivers of the microbial community (Figure 8 in the paper)
* *Figure-9-Overview_of_edge_counts_for_selected_taxa.R* - produces co-occurence networks using output from CoNet and summarizes edges within the networks (Figure 9 in the paper)
* *Sub_networks_of_the_FL_and_PA_fractions_in_the_EGC_and_WSC_regions.R* - produces sub-networks for each fraction and region (Figure 10 in the paper)

## Author:
Eduard Fadeev([eddie.fadeev-at-outlook.com](eddie.fadeev@outlook.com)) 

## Software Versions:
R version: 3.5.0\
RStudio version: 1.1.447
