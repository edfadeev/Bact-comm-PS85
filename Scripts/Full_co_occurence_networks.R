#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#load libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("cowplot"); packageVersion("cowplot")
library("dplyr"); packageVersion("dplyr")
library("stringr"); packageVersion("stringr")
library("igraph"); packageVersion("igraph")
library("parallel"); packageVersion("parallel")

#load scripts
source('./network_effect_size.R')

#load phyloseq objects
BAC_pruned <- readRDS("../Data/BAC_pruned.rds")
EUK_pruned <- readRDS("../Data/EUK_pruned.rds")

#load network tables after Co-Net analysis
FL.nodes <- read.csv("../Data/network/from_Cytoscape/DCM/FL_nodes.csv", header=T, as.is=T)
FL.links <- read.csv("../Data/network/from_Cytoscape/DCM/FL_edges.csv", header=T, as.is=T)
PA.nodes <- read.csv("../Data/network/from_Cytoscape/DCM/PA_nodes.csv", header=T, as.is=T)
PA.links <- read.csv("../Data/network/from_Cytoscape/DCM/PA_edges.csv", header=T, as.is=T)

#merged taxa of both domains
EUK_adjusted <- tax_table(EUK_pruned)[,c(1,3:6,6)]
merged_tax_table <- as.data.frame(rbind(tax_table(BAC_pruned),EUK_adjusted))
merged_tax_table$name <- gsub("_", "-", rownames(merged_tax_table))
merged_tax_table$Order <-  gsub("Incertae Sedis|_unclassified", "_u", merged_tax_table$Order) 

#####################################
#Build co-occurance Networks 
####################################
#FL
#use only edges supported by both methods
FL.links %>% filter (method_number >= 2) -> FL.links
FL.links.parsed <- as.data.frame(str_split_fixed(FL.links$Label, "->", 2))
FL.links.parsed$weight <- FL.links$qval
FL.links.parsed$type <- FL.links$interactionType
FL.links.parsed <- FL.links.parsed[FL.links.parsed$V1 %in% merged_tax_table$name,]
FL.links.parsed <- FL.links.parsed[FL.links.parsed$V2 %in% merged_tax_table$name,]
colnames(FL.links.parsed) <- c("from", "to","weight","type")

#parse taxa for network vertices
FL.nodes <- FL.nodes[which(FL.nodes$name %in% FL.links.parsed$from |
                                FL.nodes$name %in% FL.links.parsed$to),]
FL.nodes.parsed<- data.frame(name = FL.nodes$name)
FL.nodes.taxa <- plyr::join(FL.nodes.parsed, merged_tax_table, by="name", type = "left")

#make networks
FL.network <- graph_from_data_frame(d=FL.links.parsed, vertices=FL.nodes.taxa, directed=F) 

#PA
#use only edges supported by both methods
PA.links %>% filter (method_number >= 2) -> PA.links
PA.links.parsed <- as.data.frame(str_split_fixed(PA.links$Label, "->", 2))
PA.links.parsed$weight <- PA.links$qval
PA.links.parsed$type <- PA.links$interactionType
PA.links.parsed <- PA.links.parsed[PA.links.parsed$V1 %in% merged_tax_table$name,]
PA.links.parsed <- PA.links.parsed[PA.links.parsed$V2 %in% merged_tax_table$name,]
colnames(PA.links.parsed) <- c("from", "to","weight","type")

#parse taxa for network vertices
PA.nodes <- PA.nodes[which(PA.nodes$name %in% PA.links.parsed$from |
                                     PA.nodes$name %in% PA.links.parsed$to),]
PA.nodes.parsed<- data.frame(name = PA.nodes$name)
PA.nodes.taxa <- plyr::join(PA.nodes.parsed, merged_tax_table, by="name", type = "left")

#make networks
PA.network <- graph_from_data_frame(d=PA.links.parsed, vertices=PA.nodes.taxa, directed=F) 

#remove loops
FL.network <- igraph::simplify(FL.network, remove.loops = T, remove.multiple = F)
PA.network <- igraph::simplify(PA.network, remove.loops = T, remove.multiple = F)

# stack both fractions into a list
networks <- list()
networks[[1]] <- FL.network
networks[[2]] <- PA.network

#####################################
#Full network diagrams
#####################################
par(mfrow=c(1,2))
for (i in 1:2){
  #nodes colour by domain
  V(networks[[i]])$color <- V(networks[[i]])$Domain
  V(networks[[i]])$color= gsub("Bacteria","tomato",V(networks[[i]])$color)
  V(networks[[i]])$color= gsub("Eukaryota","gold",V(networks[[i]])$color)
  #Nodes names by taxonomy
  V(networks[[i]])$name <- V(networks[[i]])$Order 
  #V(networks[[i]])$name <- gsub("_unclassified","",V(networks[[i]])$name)
  
  #nodes width by degree
  V(networks[[i]])$width<- degree(networks[[i]], mode = "total")
 
  #edge colour by interraction
  E(networks[[i]])$colour <- E(networks[[i]])$type
  E(networks[[i]])$colour= gsub("copresence","green",E(networks[[i]])$colour)
  E(networks[[i]])$colour= gsub("mutualExclusion","blue",E(networks[[i]])$colour)

  #plot
  plot(networks[[i]] ,layout= layout_with_fr, 
       edge.size = (E(networks[[i]])$weight)/100, 
       edge.color= E(networks[[i]])$colour,
       vertex.size= 5+V(networks[[i]])$width,
       vertex.label =V(networks[[i]])$name,
       rescale=T,
       edge.curved=.1
       #vertex.label = NA
  )
  #add legend
  legend(x=-1.5, y=-1.1, c("Bacteria","Eukaryota"), pch=21,title= "Nodes",
         col="#777777", pt.bg=c("tomato","gold"), pt.cex=2, cex=.8, bty="n", ncol=1)
  
  legend(x=-0.7, y=-1.1, c("copresence","mutualExclusion"), pch=22, title = "Edges",
         col="#777777", pt.bg=c("green","blue"), pt.cex=2, cex=.8, bty="n", ncol=1)
  
}

#####################################
#Supplementary Table 2: Properties of the co-occurrence networks 
####################################
gs <- mclapply(networks, function(x) {
  # renaming and reversing the weight attribute 
  E(x)$wt <- max(E(x)$weight) - E(x)$weight + 1
  remove.edge.attribute(x, 'weight')
  # Centrality
  V(x)$degree      <- degree(x, mode = "total")
  V(x)$indegree    <- degree(x, mode = "in")
  V(x)$outdegree   <- degree(x, mode = "out")
  V(x)$betweenness <- betweenness(x, weights=E(x)$wt)
  V(x)$evcent      <- evcent(x)$vector
  #V(x)$closeness   <- closeness(x)
    # Local position
  V(x)$effsize     <- effective.size(x, mode = "all")
  V(x)$constraint  <- constraint(x)
  # Clustering
  com <- cluster_edge_betweenness(x,weights = E(x)$wt, modularity= TRUE)
  V(x)$memb        <- com$membership
  # Whole network
  x <- set.graph.attribute(x, "size", vcount(x))
  x <- set.graph.attribute(x, "edgecount", ecount(x))
  x <- set.graph.attribute(x, "density", edge_density(x))
  x <- set.graph.attribute(x, "mean.degree", mean(degree(x, mode = "total")))
  x <- set.graph.attribute(x, "no.clusters", length(com))
  x <- set.graph.attribute(x, "modularity", mean(com$modularity))
  x <- set.graph.attribute(x, "avg.betweeness", mean(com$edge.betweenness))
  x <- set.graph.attribute(x, "avgpathlength", average.path.length(x))
  x <- set.graph.attribute(x, "diameter", diameter(x, directed = FALSE, weights = E(x)$wt))
  x <- set.graph.attribute(x, "pos.edges", length(E(x)$type[E(x)$type == "copresence"]))
  x <- set.graph.attribute(x, "pos.edges(%)", format(round(length(E(x)$type[E(x)$type == "copresence"])/length(E(x)$type)*100,2), nsmall = 2))
  x <- set.graph.attribute(x, "neg.edges", length(E(x)$type[E(x)$type == "mutualExclusion"]))
  return(x)
})

#graph summary
gstats <- do.call('rbind', lapply(gs, function(y) {
  ga <- list.graph.attributes(y)
  sapply(ga, function(x) {
    as.numeric(get.graph.attribute(y, x))
  })
}))

network.stats <- gstats[,c("size","edgecount", "pos.edges", "pos.edges(%)", "neg.edges","no.clusters", "modularity", "diameter", "mean.degree","avgpathlength","avg.betweeness")]
network.stats <- as.table(t(network.stats))
colnames(network.stats) <- c("FL","PA")

#print table
network.stats
