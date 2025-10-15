library(dplyr)
library(RColorBrewer)
library(ggspatial)
library(ggplot2)
library(ggcorrplot)
library(cowplot)
library(ggpubr)
library(readr)
library(phyloseq)
library(patchwork)
library(purrr)
# Data prep ---------------------------------------------------------------
level <- 'asv'
#level <- 'species'
#level <- 'genus'

alpha <- 
  read_rds(paste0('./data/ati-alpha-diversity-metrics-', level, '.rds')) %>% 
  dplyr::mutate(site_id = paste0(Site, '_', Year)) %>% 
  dplyr::select(site_id, Year, Microbial_Richness, Microbial_Chao1,
                Microbial_Shannon_Diversity, Microbial_Phylogenetic_Diversity,
                Microbial_Evenness) %>% 
  dplyr::distinct(site_id, .keep_all = TRUE)

pcoa <- 
  read_rds(paste0('./data/ati-pcoa-bray-3-axes-metadata-', level, '.rds')) %>% 
  dplyr::mutate(site_id = paste0(Site, '_', Year)) %>% 
  dplyr::select(site_id, PCoA1) %>% 
  dplyr::distinct(site_id, .keep_all = TRUE)

meta <- 
  readr::read_rds(paste0('./data/ati-physeq-rrfy-', level, '.rds')) %>% 
  sample_data() %>% 
  data.frame() %>% 
  dplyr::mutate(site_id = paste0(Site, '_', Year)) %>% 
  dplyr::select(site_id, Percent_C, Percent_N, C_to_N_ratio, Phosphate,
                Silicate, Nitrite_plus_Nitrate, dN15, dC13,
                Ammonia) %>% 
  dplyr::distinct(site_id, .keep_all = TRUE)

ati <- alpha %>% 
  left_join(pcoa, join_by(site_id)) %>% 
  left_join(meta, join_by(site_id)) %>% 
  dplyr::select(-site_id) %>% 
  dplyr::rename('Species Richness' = 'Microbial_Richness',
                'Chao1 Richness' = 'Microbial_Chao1',
                'Shannon Diversity' = 'Microbial_Shannon_Diversity',
                'Phylogenetic Diversity' = 'Microbial_Phylogenetic_Diversity',
                'Evenness' = 'Microbial_Evenness',
                'PCoA 1' = 'PCoA1',
                '% C' = 'Percent_C',
                '% N' = 'Percent_N',
                'C-to-N Ratio' = 'C_to_N_ratio',
                'Nitrite+Nitrate' = 'Nitrite_plus_Nitrate',
                'delta 15N' = 'dN15',
                'delta 13C' = 'dC13') %>% 
  tidyr::drop_na()

correlations <- 
  ati %>%
  group_by(Year) %>%
  summarise(correlation_matrix = list(cor(pick(2:ncol(ati)-1))), 
            .groups = "drop")


plots <- map2(
  correlations$correlation_matrix, # List of correlation matrices
  correlations$Year,               # Corresponding years
  ~ ggcorrplot::ggcorrplot(.x, 
                           title = .y,
                           method = "square", 
                           hc.order = FALSE,
                           type = 'lower',
                           legend.title = 'Correlation',
                           outline.color = 'lightgrey',
                           colors = c("#6D9EC1", "white", "#E46726"),
                           ggtheme = theme_minimal()) # Generate plot for each year
)

combined_plot <- wrap_plots(plots, ncol = 3) + # Adjust `ncol` or `nrow` as needed
  plot_layout(guides = "collect")  # Share the legend
  
print(combined_plot)

ggsave(combined_plot, 
       file='./output/plots/corr_plots_21.22.23.png', 
       width=14.8, 
       height=5.3,
       dpi=1000)
