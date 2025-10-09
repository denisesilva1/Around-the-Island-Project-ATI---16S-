#load libraries
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(sf)
library(RColorBrewer)
library(ggspatial)

level <- 'family'; level2 <- 'Family'

physeq_rrfy <- 
  readr::read_rds(paste0('./data/ati-physeq-rrfy-', level2, '.rds')) %>% 
  #transform_sample_counts(function(x) (x / sum(x))) %>% 
  subset_taxa(Family == "Vibrionaceae" | 
                Family == "Nitrincolaceae" | 
                Family == "Pseudoalteromonadaceae" | 
                Family == "Moraxellaceae" |
                Family == "Halomonadaceae" |
                Family == "Xenococcaceae")

taxa <- tax_table(physeq_rrfy)@.Data %>% 
  data.frame() %>% 
  dplyr::mutate(asv = row.names(.)) %>% 
  dplyr::select(asv, Family)

abundances <- 
  phyloseq::sample_data(physeq_rrfy) %>% 
  as_tibble() %>% 
  dplyr::select(Site, Latitude, Longitude, Year) %>% 
  bind_cols(
    phyloseq::otu_table(physeq_rrfy) %>% 
      t()) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Verify taxa names and order
#colnames(abundances)[4:ncol(abundances)] 
colnames(abundances)[grep(paste(taxa$asv, collapse = "|"),
                     colnames(abundances))] <- taxa$Family

#Geospatial layers:
bbox <- st_bbox(abundances)
moorea <- st_read('./data/map-data/moorea.shp')

#set theme parameters
mytheme <- theme(legend.position="bottom",
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 legend.text=element_text(size=12),
                 legend.title=element_text(size=12)) 


m1 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=abundances, aes(color=Vibrionaceae+1), size=3, alpha=1) +
  theme_bw() +
  mytheme + 
  scale_color_distiller(palette='Blues', 'Vibrionaceae', direction = 1, trans = "log10") +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m1
ggsave(m1, filename='./output/maps/Vibrionaceae_21.22.23_family.png', 
       width = 11, height = 4, dpi=600)

m2 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=abundances, aes(color=Nitrincolaceae+1), size=3, alpha=1) +
  theme_bw() +
  mytheme + 
  scale_color_distiller(palette='BuGn', 'Nitrincolaceae', direction = 1, trans = "log10") +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m2
ggsave(m2, filename='./output/maps/Nitrincolaceae_21.22.23_family.png', 
       width = 11, height = 4, dpi=600)

m3 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=abundances, aes(color=Pseudoalteromonadaceae+1), size=3, alpha=1) +
  theme_bw() +
  mytheme + 
  scale_color_distiller(palette='BuPu', 'Pseudoalteromonadaceae', direction = 1, trans = "log10") +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m3
ggsave(m3, filename='./output/maps/Pseudoalteromonadaceae_21.22.23_family.png', 
       width = 11, height = 4, dpi=600)

m4 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=abundances, aes(color=Moraxellaceae+1), size=3, alpha=1) +
  theme_bw() +
  mytheme + 
  scale_color_distiller(palette='GnBu', 'Moraxellaceae', direction = 1, trans = "log10") +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m4
ggsave(m4, filename='./output/maps/Moraxellaceae_21.22.23_family.png', 
       width = 11, height = 4, dpi=600)

m5 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=abundances, aes(color=Halomonadaceae+1), size=3, alpha=1) +
  theme_bw() +
  mytheme + 
  scale_color_distiller(palette='Oranges', 'Halomonadaceae', direction = 1, trans = "log10") +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m5
ggsave(m5, filename='./output/maps/Halomonadaceae_21.22.23_familiy.png', 
       width = 11, height = 4, dpi=600)

m6 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=abundances, aes(color=Xenococcaceae+1), size=2, alpha=1) +
  theme_bw() +
  mytheme + 
  scale_color_distiller(palette='Oranges', 'Xenococcaceae', direction = 1, trans = "log10") +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m6
ggsave(m6, filename='./output/maps/Xenococcaceae_21.22.23_familiy.png', 
       width = 11, height = 4, dpi=600)



