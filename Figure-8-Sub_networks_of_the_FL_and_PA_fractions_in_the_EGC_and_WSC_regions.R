#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#load libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("cowplot"); packageVersion("cowplot")
library("dplyr"); packageVersion("dplyr")
library("igraph"); packageVersion("igraph")
library("parallel"); packageVersion("parallel")

#rerun enrichment test
source("./Scripts/Enriched_bacterial_daOTU.R")
source("./Scripts/Full_co_occurence_networks.R")

#####################################
#Figure 7: Sub-networks of taxon-taxon associations
#####################################
##prepare nodes for networks
#FL
#assign the enrichment for each bacterial taxa
FL.nodes.taxa.enr <- plyr::join(FL.nodes.taxa, BAC_daOTU, by="name", type = "left")
FL.nodes.taxa.enr$region <- as.character(FL.nodes.taxa.enr$region)
FL.nodes.taxa.enr[which(FL.nodes.taxa.enr$Domain=="Bacteria" & is.na(FL.nodes.taxa.enr$region)),c("region")] <- "none"
FL.nodes.taxa.enr[which(FL.nodes.taxa.enr$Domain=="Eukaryota"),c("region")] <- "UNDEF"

#PA
#assign the enrichment for each bacterial taxa
PA.nodes.taxa.enr <- plyr::join(PA.nodes.taxa, BAC_daOTU, by="name", type = "left")
PA.nodes.taxa.enr$region <- as.character(PA.nodes.taxa.enr$region)
PA.nodes.taxa.enr[which(PA.nodes.taxa.enr$Domain=="Bacteria" & is.na(PA.nodes.taxa.enr$region)),c("region")] <- "none"
PA.nodes.taxa.enr[which(PA.nodes.taxa.enr$Domain=="Eukaryota"),c("region")] <- "UNDEF"


##prepare edges with taxonomy for networks
#assign taxa to edges for further networks
merged_taxa_from <- merged_tax_table
colnames(merged_taxa_from)[7] <- c("from")
merged_taxa_to <- merged_tax_table
colnames(merged_taxa_to)[7] <- c("to")

##loop the network building
fraction_FL <- c("FL.EGC", "FL.WSC", "none")
fraction_PA <- c("PA.EGC", "PA.WSC", "none")
nodes_parsed <- list()
links_parsed <- list()
nets_FL <- list()
nets_PA<- list()
edges_FL <- list()
nodes_FL <- list()
edges_PA <- list()
nodes_PA <- list()

#FL
for (i in 1:3){
  nodes_parsed[[i]] <- as.data.frame(subset(FL.nodes.taxa.enr, region == as.character(fraction_FL[i]) | region == "UNDEF"))
  links_parsed[[i]]  <- as.data.frame(FL.links.parsed[which(FL.links.parsed$from %in% nodes_parsed[[i]]$name & FL.links.parsed$to %in% nodes_parsed[[i]]$name),])
  
  FL_from <- left_join(links_parsed[[i]], merged_taxa_from, by="from", type = "left")
  FL_to <- left_join(links_parsed[[i]], merged_taxa_to, by="to", type = "left")
  
  FL_merged_edges <- data.frame(FL_from[,c("Order")])
  FL_merged_edges <- cbind(FL_merged_edges,FL_to[,c("Order","weight","type")])
  colnames(FL_merged_edges) <- c("from","to","weight","type")
  
  FL_merged_edges_pos_agg <- aggregate(weight~from+to+type, FL_merged_edges[FL_merged_edges$type=="copresence",],FUN = length)
  
  nodes_parsed_agg <-aggregate(name~Domain+Phylum+Class+Order, nodes_parsed[[i]],FUN = length) 
  
  colnames(nodes_parsed_agg) <- c("Domain", "Phylum","Class","name","OTUs")
  nodes_parsed_agg <- nodes_parsed_agg[,c("name", "OTUs", "Domain")]
  
  nodesFL_parsed_taxa_enr_FL_agg <- nodes_parsed_agg[which(nodes_parsed_agg$name %in% FL_merged_edges_pos_agg$from |
                                                             nodes_parsed_agg$name %in% FL_merged_edges_pos_agg$to),]
  
  netFL_pos <- graph_from_data_frame(d=FL_merged_edges_pos_agg, 
                                     vertices=nodesFL_parsed_taxa_enr_FL_agg, directed=F) 
  
  nets_FL[[i]] <-igraph::simplify(netFL_pos, remove.loops = T, remove.multiple = T)
  
  V(nets_FL[[i]])$degree <- igraph::degree(nets_FL[[i]], mode="all", normalized = FALSE)
  V(nets_FL[[i]])$strength <- igraph::strength(nets_FL[[i]], mode="all")
  #nodes colour
  V(nets_FL[[i]])$color <- V(nets_FL[[i]])$Domain
  V(nets_FL[[i]])$color= gsub("Bacteria","tomato",V(nets_FL[[i]])$color)
  V(nets_FL[[i]])$color= gsub("Eukaryota","gold",V(nets_FL[[i]])$color)
  
  #edges color
  #E(nets_FL[[i]])$color <- ifelse(E(nets_FL[[i]])$color > 0, "green", "gray") 
  
  edges_FL[[i]] <- as_data_frame(nets_FL[[i]], what = c("edge"))
  nodes_FL[[i]] <- as_data_frame(nets_FL[[i]], what = c("vertices"))
}

#PA
for (i in 1:3){
  nodes_parsed[[i]] <- as.data.frame(subset(PA.nodes.taxa.enr, region == as.character(fraction_PA[i]) | region == "UNDEF"))
  links_parsed[[i]]  <- as.data.frame(PA.links.parsed[which(PA.links.parsed$from %in% nodes_parsed[[i]]$name & PA.links.parsed$to %in% nodes_parsed[[i]]$name),])
  
  PA_from <- left_join(links_parsed[[i]], merged_taxa_from, by="from", type = "left")
  PA_to <- left_join(links_parsed[[i]], merged_taxa_to, by="to", type = "left")
  
  PA_merged_edges <- data.frame(PA_from[,c("Order")])
  PA_merged_edges <- cbind(PA_merged_edges,PA_to[,c("Order","weight","type")])
  colnames(PA_merged_edges) <- c("from","to","weight","type")
  
  PA_merged_edges_pos_agg <- aggregate(weight~from+to+type, PA_merged_edges[PA_merged_edges$type=="copresence",],FUN = length)
  
  nodes_parsed_agg <-aggregate(name~Domain+Phylum+Class+Order, nodes_parsed[[i]],FUN = length) 
  
  colnames(nodes_parsed_agg) <- c("Domain", "Phylum","Class","name","OTUs")
  nodes_parsed_agg <- nodes_parsed_agg[,c("name", "OTUs", "Domain")]
  
  nodesPA_parsed_taxa_enr_PA_agg <- nodes_parsed_agg[which(nodes_parsed_agg$name %in% PA_merged_edges_pos_agg$from |
                                                             nodes_parsed_agg$name %in% PA_merged_edges_pos_agg$to),]
  
  netPA_pos <- graph_from_data_frame(d=PA_merged_edges_pos_agg, 
                                     vertices=nodesPA_parsed_taxa_enr_PA_agg, directed=F) 
  
  nets_PA[[i]] <-igraph::simplify(netPA_pos, remove.loops = T, remove.multiple = T)
  
  V(nets_PA[[i]])$degree <- igraph::degree(nets_PA[[i]], mode="all", normalized = FALSE)
  V(nets_PA[[i]])$strength <- igraph::strength(nets_PA[[i]], mode="all")
  #nodes colour
  V(nets_PA[[i]])$color <- V(nets_PA[[i]])$Domain
  V(nets_PA[[i]])$color= gsub("Bacteria","tomato",V(nets_PA[[i]])$color)
  V(nets_PA[[i]])$color= gsub("Eukaryota","gold",V(nets_PA[[i]])$color)
  
  #edges color
  #E(nets_PA[[i]])$color <- ifelse(E(nets_PA[[i]])$color > 0, "green", "gray") 
  
  edges_PA[[i]] <- as_data_frame(nets_PA[[i]], what = c("edge"))
  nodes_PA[[i]] <- as_data_frame(nets_PA[[i]], what = c("vertices"))
}


#Plot the networks
png('./Figures/Figure-8-sub-networks.png',width = 30, height = 30, res=300,units = "cm")
par(mfrow=c(2,2))
for (i in 1:2){
  l <- layout_in_circle(nets_FL[[i]])
  #l <- norm_coords(l, ymin=-1.1, ymax=1.1, xmin=-1.1, xmax=1.1)
  plot(nets_FL[[i]], layout = l, 
       vertex.size = 6+V(nets_FL[[i]])$OTUs/3, 
       #vertex.size = V(nets_FL[[i]])$strength*0.3, 
       vertex.size = V(nets_FL[[i]])$degree, 
       edge.width= 2+(E(nets_FL[[i]])$weight)/3,
       vertex.label =V(nets_FL[[i]])$name,
       rescale=T
       #vertex.label = NA
  )
}

for (i in 1:2){
  l <- layout_in_circle(nets_PA[[i]])
  #l <- norm_coords(l, ymin=-1.1, ymax=1.1, xmin=-1.1, xmax=1.1)
  plot(nets_PA[[i]], layout = l, 
       vertex.size =6+V(nets_PA[[i]])$OTUs/3,
       #vertex.size = V(nets_PA[[i]])$strength*0.3,
       #vertex.size = V(nets_PA[[i]])$degree,
       edge.width= 2+(E(nets_PA[[i]])$weight)/3,
       vertex.label =V(nets_PA[[i]])$name,
       rescale=T
       #vertex.label = NA
  )
}

dev.off()

#####################################
#remove all temporary datasets 
####################################
detach("package:phyloseq", unload=TRUE)
detach("package:parallel", unload=TRUE)
detach("package:ggplot2", unload=TRUE)
detach("package:dplyr", unload=TRUE)
detach("package:cowplot", unload=TRUE)
detach("package:igraph", unload=TRUE)

rm(list=ls(all=TRUE))

sessionInfo()



