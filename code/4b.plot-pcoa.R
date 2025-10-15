library(ggplot2)
library(dplyr)
library(spmodel)
library(sf)
library(car)
library(vegan)

shore_colors <- c('#1f78b4', '#fb9a99', '#33a02c')
habitat_colors <- brewer.pal(5, "Set1")

level <- 'asv'

dis <- 'bray'; dis2 = "Bray-Curtis"
#dis <- 'unifrac'; dis2 = "Weighted Unifrac"

physeq <- 
  readr::read_rds(paste0('./data/ati-physeq-nonrrfy-', level, '.rds'))

pcoa <- 
  readr::read_rds(paste0('./data/ati-pcoa-', dis, '-object-', level, '.rds')) 

plot_ordination(physeq, 
                pcoa, 
                type = "Site", 
                color = "Habitat", 
                title = paste0(dis2, " PCoA")) +
  geom_point(size = 4, aes(shape = Island_shore)) +
  scale_color_manual(values=habitat_colors) +
  facet_wrap(~Year) +
  labs(shape = "Island Shore") +
  theme_bw() 

ggplot2::ggsave(file = paste0('./output/plots/pcoa-plot-', 
                              dis, '-', 
                              level, '.png'),
                dpi = 600,
                width = 300,
                units = 'mm')



