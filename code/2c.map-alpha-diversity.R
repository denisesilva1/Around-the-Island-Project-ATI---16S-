library(sf)
library(dplyr)
library(RColorBrewer)
library(ggspatial)
library(ggplot2)

#Prep colors
# shore_colors <- c('#1f78b4', '#fb9a99', '#33a02c')
# habitat_colors <- brewer.pal(5, "Set1")

level <- 'asv'

alpha_div <- 
  readr::read_rds(
    paste0('./data/ati-alpha-diversity-metrics-', level, '.rds'))

#Geospatial layers
site_data <- st_as_sf(alpha_div, coords = c("Longitude", "Latitude"))
bbox <- st_bbox(site_data)
moorea <- st_read('./data/map-data/moorea.shp')
st_crs(site_data) <- st_crs(moorea)

#set theme parameters
mytheme <- theme(legend.position="bottom",
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 legend.text=element_text(size=12),
                 legend.title=element_text(size=12)) 


m2 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data, aes(color=Microbial_Richness), size=3, alpha=1) +
  theme_bw() +
  mytheme + 
  scale_color_distiller(palette='Blues', 'Microbial Richness', direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m2
ggsave(m2, 
       filename = paste0('./output/maps/microbial_richness_',
                         level,
                         '_21.22.23.png'), 
       width = 11, height = 4, dpi=600)

m3 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data, aes(color=Microbial_Shannon_Diversity), size=3, alpha=1) +
  theme_bw() +
  mytheme + scale_color_distiller(palette='BuGn', 'Microbial Shannon Diversity', direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m3
ggsave(m3, 
       filename = paste0('./output/maps/microbial_shannon_',
                         level,
                         '_21.22.23.png'), 
       width = 11, height = 4, dpi=600)

m4 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data, aes(color=Microbial_Phylogenetic_Diversity), size=3, alpha=1) +
  theme_bw() +
  mytheme + 
  mytheme + scale_color_distiller(palette='BuPu', 'Microbial Phylogenetic Diversity', direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m4
ggsave(m4, 
       filename = paste0('./output/maps/microbial_phylo_diversity_',
                         level,
                         '_21.22.23.png'), 
       width = 11, height = 4, dpi=600)

m5 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data, aes(color=Microbial_Evenness), size=3, alpha=1) +
  theme_bw() +
  mytheme + 
  mytheme + scale_color_distiller(palette='GnBu', 'Microbial Evenness', direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m5
ggsave(m5, 
       filename = paste0('./output/maps/microbial_evennes_',
                         level,
                         '_21.22.23.png'), 
       width = 11, height = 4, dpi=600)


m6 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data, aes(color=Microbial_Evenness), size=3, alpha=1) +
  theme_bw() +
  mytheme + 
  mytheme + scale_color_distiller(palette='Oranges', 'Microbial Chao1', direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m6
ggsave(m6, 
       filename = paste0('./output/maps/microbial_chao1_',
                         level,
                         '_21.22.23.png'), 
       width = 11, height = 4, dpi=600)

