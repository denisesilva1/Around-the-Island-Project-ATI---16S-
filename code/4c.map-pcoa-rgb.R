library(ggplot2)
library(dplyr)
library(spmodel)
library(sf)
library(car)
library(vegan)

pcoa.colors <- function(pcoa){
  min0 <- min(pcoa)
  rn0 <- diff(range(pcoa))
  scale.to.color <- function(x, min0, rn0){
    format.hexmode(round((x-min0)/rn0*255))
  }
  sc.out <- apply(pcoa, 2, scale.to.color, 
                  min0 = min0, rn0 = rn0)
  col0 <- apply(sc.out, 1, 
                function(x) paste("#", 
                                  paste(x, collapse = ""),
                                  sep = ""))
  
  return(col0)
}

#Geospatial layers:
bbox <- st_bbox(ati.pcoa)
moorea <- st_read('./data/map-data/moorea.shp')

#set theme parameters
mytheme <- theme(legend.position="bottom",
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 legend.text=element_text(size=12),
                 legend.title=element_text(size=12)) 


# Create RGB from 3 PCoA axes ---------------------------------------------
level <- 'asv'
#level <-'species'
#level <- 'genus'

dis <- 'bray'
#dis <- 'unifrac'

ati.pcoa <- 
  readr::read_rds(paste0('./data/ati-pcoa-', dis, '-3-axes-metadata-', 
                         level, 
                         '.rds')) %>% 
  mutate(Year = as.factor(Year)) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs=4326)

pcoa.cols <- ati.pcoa %>% 
  st_drop_geometry()  %>% 
  dplyr::select(c(PCoA1, PCoA2, PCoA3))

pcoa.hex <- pcoa.colors(pcoa.cols)

ati.pcoa <- cbind(ati.pcoa, pcoa.hex)

p <- 
  ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=ati.pcoa, 
          aes(color=pcoa.hex), 
          size=2, show.legend = F, alpha=0.65) +
  scale_colour_identity() +
  theme_bw() + mytheme + 
  #guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)

ggsave(filename = paste0('./output/maps/PCoA_', dis, '_RGB_21.22.23_', level, '.png'), 
       width = 11, height = 4, dpi=600)

readr::write_rds(p, 
                 file = './data/combined_figures/combined_pcoa_rgb_map_by_year.rds')

library(ggplot2)
library(dplyr)
library(spmodel)
library(sf)
library(car)
library(vegan)

pcoa.colors <- function(pcoa){
  min0 <- min(pcoa)
  rn0 <- diff(range(pcoa))
  scale.to.color <- function(x, min0, rn0){
    format.hexmode(round((x-min0)/rn0*255))
  }
  sc.out <- apply(pcoa, 2, scale.to.color, 
                  min0 = min0, rn0 = rn0)
  col0 <- apply(sc.out, 1, 
                function(x) paste("#", 
                                  paste(x, collapse = ""),
                                  sep = ""))
  
  return(col0)
}

#Geospatial layers:
bbox <- st_bbox(ati.pcoa)
moorea <- st_read('./data/map-data/moorea.shp')

#set theme parameters
mytheme <- theme(legend.position="bottom",
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 legend.text=element_text(size=12),
                 legend.title=element_text(size=12)) 


# Create RGB from 3 PCoA axes ---------------------------------------------
level <- 'asv'
#level <-'species'
#level <- 'genus'

dis <- 'bray'
#dis <- 'unifrac'

ati.pcoa <- 
  readr::read_rds(paste0('./data/ati-pcoa-', dis, '-3-axes-metadata-', 
                         level, 
                         '.rds')) %>% 
  mutate(Year = as.factor(Year)) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs=4326)

pcoa.cols <- ati.pcoa %>% 
  st_drop_geometry()  %>% 
  dplyr::select(c(PCoA1, PCoA2, PCoA3))

pcoa.hex <- pcoa.colors(pcoa.cols)

ati.pcoa <- cbind(ati.pcoa, pcoa.hex)

p <- 
  ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=ati.pcoa, 
          aes(color=pcoa.hex), 
          size=2, show.legend = F, alpha=0.65) +
  scale_colour_identity() +
  theme_bw() + mytheme + 
  #guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)

p
