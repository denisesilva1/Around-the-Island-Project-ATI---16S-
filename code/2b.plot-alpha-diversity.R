library(sf)
library(dplyr)
library(RColorBrewer)
library(ggspatial)
library(ggplot2)
library(ggbeeswarm)

#Prep colors
shore_colors <- c('#1f78b4', '#fb9a99', '#33a02c')
habitat_colors <- brewer.pal(5, "Set1")

level <- 'asv'

alpha_div <- 
  readr::read_rds(
    paste0('./data/ati-alpha-diversity-metrics-', level, '.rds')) 

alpha_div$Habitat <- recode(
  alpha_div$Habitat,
  "Fringing reef" = "Fringing\nreef",
  "Mid lagoon" = "Mid\nlagoon",
  "Reef crest" = "Reef\ncrest",
  "Reef pass" = "Reef\npass")

# richness ----------------------------------------------------------------
ggplot(alpha_div, aes(x = Habitat, 
                      y = Microbial_Richness, 
                      color = Habitat)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.9, size=0.8) +  # Boxplot overlaid
  geom_quasirandom(size = 2.5, alpha = 0.3, stroke = NA) +
  #geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.5) +  # Boxplot overlaid
  scale_fill_manual(values = habitat_colors) +
  scale_color_manual(values = habitat_colors) +
  labs(x = "Habitat Type", y = "Microbial Richness") +
  facet_wrap(~Year) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file = paste0('./output/plots/richness-by-habitat-boxplots-', level, '.png'),
       units = 'mm',
       width = 230,
       height = 85,
       dpi = 600)

# chao1 ----------------------------------------------------------------
ggplot(alpha_div, aes(x = Habitat, 
                      y = Microbial_Chao1, 
                      color = Habitat)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.9, size=0.8) +  # Boxplot overlaid
  geom_quasirandom(size = 2.5, alpha = 0.3, stroke = NA) +
  #geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.5) +  # Boxplot overlaid
  scale_fill_manual(values = habitat_colors) +
  scale_color_manual(values = habitat_colors) +
  labs(x = "Habitat Type", y = "Chao1 Richness") +
  facet_wrap(~Year) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file = paste0('./output/plots/chao1-by-habitat-boxplots-', level, '.png'),
       units = 'mm',
       width = 230,
       height = 85,
       dpi = 600)

# Shannon ----------------------------------------------------------------
ggplot(alpha_div, aes(x = Habitat, 
                      y = Microbial_Shannon_Diversity, 
                      color = Habitat)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.9, size=0.8) +  # Boxplot overlaid
  geom_quasirandom(size = 2.5, alpha = 0.3, stroke = NA) +
  #geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.5) +  # Boxplot overlaid
  scale_fill_manual(values = habitat_colors) +
  scale_color_manual(values = habitat_colors) +
  labs(x = "Habitat Type", y = "Shannon Diversity") +
  facet_wrap(~Year) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file = paste0('./output/plots/shannon-by-habitat-boxplots-', level, '.png'),
       units = 'mm',
       width = 230,
       height = 85,
       dpi = 600)

# Phylogenetic  ----------------------------------------------------------------
ggplot(alpha_div, aes(x = Habitat, 
                      y = Microbial_Phylogenetic_Diversity, 
                      color = Habitat)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.9, size=0.8) +  # Boxplot overlaid
  geom_quasirandom(size = 2.5, alpha = 0.3, stroke = NA) +
  #geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.5) +  # Boxplot overlaid
  scale_fill_manual(values = habitat_colors) +
  scale_color_manual(values = habitat_colors) +
  labs(x = "Habitat Type", y = "Phylogenetic Diversity") +
  facet_wrap(~Year) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file = paste0('./output/plots/phylogenetic-by-habitat-boxplots-', level, '.png'),
       units = 'mm',
       width = 230,
       height = 85,
       dpi = 600)

# Evenness----------------------------------------------------------------
ggplot(alpha_div, aes(x = Habitat, 
                      y = Microbial_Evenness, 
                      color = Habitat)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.9, size=0.8) +  # Boxplot overlaid
  geom_quasirandom(size = 2.5, alpha = 0.3, stroke = NA) +
  #geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.5) +  # Boxplot overlaid
  scale_fill_manual(values = habitat_colors) +
  scale_color_manual(values = habitat_colors) +
  labs(x = "Habitat Type", y = "Microbial Evenness") +
  facet_wrap(~Year) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file = paste0('./output/plots/evenness-by-habitat-boxplots-', level, '.png'),
       units = 'mm',
       width = 230,
       height = 85,
       dpi = 600)
