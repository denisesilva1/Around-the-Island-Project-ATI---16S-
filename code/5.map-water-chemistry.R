library(sf)
library(dplyr)
library(RColorBrewer)
library(ggspatial)
library(ggplot2)

#Prep colors
shore_colors <- c('#1f78b4', '#fb9a99', '#33a02c')
habitat_colors <- brewer.pal(5, "Set1")

level <- 'asv'

alpha_div <- 
  readr::read_rds(file = paste0('./data/ati-alpha-diversity-metrics-', level, '.rds'))

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

m6 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data %>% 
            filter(Nitrite_plus_Nitrate < 2.5), 
          aes(color=Nitrite_plus_Nitrate), size=3, alpha=1) +
  theme_bw() +
  #mytheme + 
  mytheme + scale_color_distiller(palette='Greens', 'Nitrate+Nitrite', direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m6
ggsave(m6, filename='./output/maps/Nitrite_plus_Nitrate_21.22.23_remove_outliers.png', 
       width = 11, height = 4, dpi=600)
readr::write_rds(m6, './data/combined_figures/Nitrite_plus_Nitrate.rds')

m7 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data %>% 
            filter(!is.na(Silicate)), 
          aes(color=Silicate), size=3, alpha=1) +
  theme_bw() +
  mytheme + scale_color_distiller(palette='Oranges', 'Silicate', direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m7
ggsave(m7, filename='./output/maps/Silicate_21.22.23.png', 
       width = 11, height = 4, dpi=600)

m8 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data %>% 
            filter(!is.na(Phosphate)), 
          aes(color=Phosphate), size=3, alpha=1) +
  theme_bw() +
  mytheme + scale_color_distiller(palette='PuRd', 'Phosphate', direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m8
ggsave(m8, filename='./output/maps/Phosphate_21.22.23.png', 
       width = 11, height = 4, dpi=600)
readr::write_rds(m8, './data/combined_figures/Phosphate.rds')


m9 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data %>% 
            filter(!is.na(C_to_N_ratio)), 
          aes(color=C_to_N_ratio), size=3, alpha=1) +
  theme_bw() +
  mytheme + scale_color_distiller(palette='Purples', 'C-to-N Ratio', direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m9
ggsave(m9, filename='./output/maps/C_to_N_ratio_21.22.23.png', 
       width = 11, height = 4, dpi=600)

m10 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data %>% 
            filter(!is.na(Ammonia)), 
          aes(color=Ammonia), size=3, alpha=1) +
  theme_bw() +
  mytheme + scale_color_distiller(palette='Reds', 'Ammonia', direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m10
ggsave(m10, filename='./output/maps/Ammonia_21.22.23.png', 
       width = 11, height = 4, dpi=600)
readr::write_rds(m10, './data/combined_figures/Ammonia.rds')

m11 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data %>% 
            filter(!is.na(Percent_N)), 
          aes(color=Percent_N), size=3, alpha=1) +
  theme_bw() +
  mytheme + scale_color_distiller(palette='YlGnBu', '% N', direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m11
ggsave(m11, filename='./output/maps/Percent_N_21.22.23.png', 
       width = 11, height = 4, dpi=600)

m12 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data %>% 
            filter(!is.na(Percent_C)), aes(color=Percent_C), size=3, alpha=1) +
  theme_bw() +
  mytheme + scale_color_distiller(palette='Greys', '% C', direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m12
ggsave(m12, filename='./output/maps/Percent_C_21.22.23.png', 
       width = 11, height = 4, dpi=600)

m13 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data %>% 
            filter(!is.na(dN15)), 
          aes(color=dN15), size=3, alpha=1) +
  theme_bw() +
  mytheme + scale_color_distiller(palette='YlGnBu', 
                                  expression(delta^15*N), 
                                  direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m13

ggsave(m13, filename='./output/maps/delta15N_21.22.23.png', 
       width = 11, height = 4, dpi=600)

m14 <-ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=site_data %>% 
            filter(!is.na(dC13)), 
          aes(color=dC13), size=3, alpha=1) +
  theme_bw() +
  mytheme + scale_color_distiller(palette='YlGnBu', 
                                  expression(delta^13*C), 
                                  direction = 1) +
  guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)
m14

ggsave(m14, filename='./output/maps/delta13C_21.22.23.png', 
       width = 11, height = 4, dpi=600)
