#load library
library(RColorBrewer)


#####################################
#Color palettes for plots
#####################################
#large colours range
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

tol21rainbow<- c("#771155", 
                 "#AA4488", 
                 "#CC99BB", 
                 "#114477", 
                 "#4477AA", 
                 "#77AADD", 
                 "#117777", 
                 "#44AAAA", 
                 "#77CCCC", 
                 "#117744", 
                 "#44AA77", 
                 "#88CCAA", 
                 "#777711", 
                 "#AAAA44", 
                 "#DDDD77", 
                 "#774411", 
                 "#AA7744", 
                 "#DDAA77", 
                 "#771122", 
                 "#AA4455", 
                 "#DD7788",
                 "#56B4E9")



phyla.col <- c("Alphaproteobacteria"= "#771155" ,
               "Acidimicrobiia"="#AA4488",
               "Betaproteobacteria"="#CC99BB", 
               "BD2-11 terrestrial group"= "#114477", 
               "Deltaproteobacteria"= "#4477AA", 
               "Flavobacteriia"= "#77AADD", 
               "Gammaproteobacteria"="#117777", 
               "Marinimicrobia (SAR406 clade)_unc"= "#34ABAA", 
               "Nitrospinia" = "#77CCCC",  
               "OM190"="#117744", 
               "Opitutae"= "#44AA77", 
               "Phycisphaerae"= "#88CCAA", 
               "Planctomycetacia"= "#777711", 
               "SAR202 clade"= "#AAAA44", 
               "Subgroup 6"= "#000000",
               "Sphingobacteriia" ="#DDDD77", 
               "SPOTSOCT00m83" = "#774411", 
               "Verrucomicrobiae" = "#AA7744")


reg_colours <-c("EGC"= "blue","WSC"= "red")

theme_plot <- theme(axis.title.x = element_text(size = 15), 
                    axis.title.y = element_text(size = 15),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 12),
                    strip.text = element_text(size = 15),
                    legend.text= element_text(size = 12),
                    plot.title = element_text(size=18, lineheight=.8, face="bold", hjust = 0.5),
                    panel.background = element_rect(fill = 'white', colour = 'white'),
                    panel.border = element_rect(color = "black", fill = NA, size = 1), 
                    strip.background = element_rect(color = "black", size = 1),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank())


gamma_col <- c("Alteromonadales" = "#000000", 
               "Arenicellales"="#E69F00", 
               "BD7-8 marine group" = "#56B4E9", 
               "Cellvibrionales"= "#009E73", 
               "E01-9C-26 marine group"= "#F0E442", 
               "Gammaproteobacteria_unc" = "#0072B2", 
               "KI89A clade" = "#D55E00", 
               "Legionellales" = "#CC79A7",
               "Oceanospirillales"= "#771155" ,
               "Pseudomonadales"="#AA4488",
               "Salinisphaerales"="#CC99BB", 
               "Thiotrichales"= "#114477", 
               "Xanthomonadales"= "#4477AA", 
               "Vibrionales"= "#77AADD")


