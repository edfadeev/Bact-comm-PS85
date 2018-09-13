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

#set plots theme
theme_set(theme_classic())

#build full networks
source("./Scripts/Full_co_occurence_networks.R")

#####################################
#Figure 7: Overview of edge counts for selected taxonomic groups in each network
####################################
vert.stats <- do.call('rbind', lapply(1:length(gs), function(x) {
  o <- get.data.frame(gs[[x]], what = 'vertices')
  o$frac <- x #1=FL, 2=PA
  return(o)
}))

edge.stats <- do.call('rbind', lapply(1:length(gs), function(x) {
  o <- get.data.frame(gs[[x]], what = 'edges')
  o$frac <- x #1=FL, 2=PA
  return(o)
}))

#separate fractions and count edges
#FL
FL.edge.stats <- edge.stats[edge.stats$frac == 1,]
FL.edge.stats_t <- FL.edge.stats[,c("to", "type")]
colnames(FL.edge.stats_t)[1] <- "from"
FL_edges.long <- rbind(FL.edge.stats[,c("from", "type")], FL.edge.stats_t)
FL_edges.long$frac <- 1

#PA
PA.edge.stats <- edge.stats[edge.stats$frac == 2,]
PA.edge.stats_t <- PA.edge.stats[,c("to", "type")]
colnames(PA.edge.stats_t)[1] <- "from"
PA_edges.long <- rbind(PA.edge.stats[,c("from", "type")], PA.edge.stats_t)
PA_edges.long$frac <- 2

#merge edges stats
edges.stats.long <- rbind(FL_edges.long,PA_edges.long)

#merge edges with nodes to one table
colnames(vert.stats)[1] <- "from"
networks.edges.stats <- merge(edges.stats.long, vert.stats, by=c("from","frac"))

#count number of edges for each node
edges.sum.by.taxa <- as.data.frame(as.list(aggregate(degree~Domain+Phylum+Class+Order+type+frac, networks.edges.stats, FUN = length )))

#filter taxa with 3 or more edges
edges.sum.lim <- as.data.frame(as.list(aggregate(degree~Domain+Phylum+Class+Order+frac, networks.edges.stats, FUN = length )))
edges.sum.lim <- edges.sum.lim[edges.sum.lim$degree>10,]
edges.sum.by.taxa <- edges.sum.by.taxa[edges.sum.by.taxa$Order %in% edges.sum.lim$Order,]


BAC_edges.p <- ggplot(edges.sum.by.taxa[edges.sum.by.taxa$Domain == "Bacteria",], 
                      aes(x=Order, y = degree, fill = type))+
                      geom_col()+ 
                      facet_grid(frac~.)+ 
                      ylab("Associations count")+#ylim(0,45)+
                      theme(legend.position = "none", axis.text.x = element_text(angle = 90))+
                      scale_fill_manual(values = c("copresence" = "#008837","mutualExclusion"= "#7b3294"))

EUK_edges.p <-ggplot(edges.sum.by.taxa[edges.sum.by.taxa$Domain == "Eukaryota",], 
                     aes(x=Order, y = degree, fill = type))+
                      geom_col()+ 
                      facet_grid(frac~.)+ 
                      ylab("Associations count")+#ylim(0,45)+
                      theme(legend.position = "none", axis.text.x = element_text(angle = 90))+
                      scale_fill_manual(values = c("copresence" = "#008837","mutualExclusion"= "#7b3294"))

png('./Figures/Figure-7-edges_counts.png',width = 30, height = 30, res=300, units = "cm")
plot_grid(BAC_edges.p, EUK_edges.p, labels = c("A", "B"), ncol = 2, align = "hv")
dev.off()

#####################################
#remove all temporary datasets 
####################################
detach("package:phyloseq", unload=TRUE)
detach("package:igraph", unload=TRUE)
detach("package:ggplot2", unload=TRUE)
detach("package:dplyr", unload=TRUE)
detach("package:cowplot", unload=TRUE)
detach("package:parallel", unload=TRUE)

rm(list=ls(all=TRUE))

sessionInfo()

