library(sf)
library(tidyverse)
library(RColorBrewer)
library(ggspatial)
library(ggplot2)
library(ggcorrplot)
library(cowplot)
library(ggtern)

# Data prep ---------------------------------------------------------------
level <- 'asv'
#level <- 'species'
#level <- 'genus'

#Prep colors
shore_colors <- c('#1f78b4', '#fb9a99', '#33a02c')
habitat_colors <- brewer.pal(5, "Set1")

vars <- c('Site', 'Year', 'Habitat', 'Phosphate', 'Silicate', 
          'Nitrite_plus_Nitrate','Ammonia', 'Microbial_Richness',
          'Microbial_Shannon_Diversity', 'Microbial_Chao1',
          'Microbial_Phylogenetic_Diversity', 'Microbial_Evenness',
          'Percent_N', 'Percent_C', 'C_to_N_ratio', 'dN15', 'dC13')

lookup <- c(`Nitrite plus Nitrate` = "Nitrite_plus_Nitrate",
            `% C` = "Percent_C",
            `% N` = "Percent_N",
            `C-to-N Ratio` = "C_to_N_ratio",
            `Microbial Chao1` = "Microbial_Chao1",
            `Microbial Richness` = "Microbial_Richness",
            `Microbial Shannon Diversity` = "Microbial_Shannon_Diversity",
            `Microbial Phylogenetic Diversity` = "Microbial_Phylogenetic_Diversity",
            `Microbial Evenness` = "Microbial_Evenness")

alpha_div <-
  readr::read_rds(file = paste0('./data/ati-alpha-diversity-metrics-', level, '.rds')) %>% 
  dplyr::select(all_of(vars)) %>% 
  dplyr::rename(all_of(lookup))

vars2 <- names(alpha_div)[-1:-3]
vars <- vars[-1:-3]

for(i in 1:length(vars2)){

  df_wide <- alpha_div %>%
    dplyr::select(Site, Year, vars2[i]) %>%
    dplyr::rename(x = vars2[i]) %>% 
    na.omit() %>% 
    dplyr::mutate(x_scaled = (x - min(x)) / (max(x) - min(x))*100) %>% 
    dplyr::select(Site, Year, x_scaled) %>%
    pivot_wider(
      names_from = Year,
      values_from = x_scaled,
      values_fn = mean) %>% 
    dplyr::right_join(alpha_div %>% 
                        dplyr::select(Site, Habitat) %>% 
                        dplyr::distinct(),
                      join_by(Site)) %>% 
    na.omit()
  
  
  m <- ggtern(
    data = df_wide,
    aes(x = `2022`, y = `2021`, z = `2023`, color = Habitat)) +
    labs(x = "2022", y = "2021", z = "2023") +
    scale_colour_manual(name = "Habitat", values = habitat_colors) +
    geom_point(alpha = 0.5, size = 2) +
    theme_linedraw() +
    theme(
      # --- Triangle outline & ticks ---
      tern.axis.line        = element_line(linewidth = 1.5),  # outer triangle (valid)
      axis.ticks            = element_line(linewidth = 1.2),
      
      # --- Gridlines ---
      tern.panel.grid.major = element_line(linewidth = 0.8),
      tern.panel.grid.minor = element_line(linewidth = 0.4),
      
      # --- Legend under plot, tightened spacing ---
      legend.position       = "bottom",
      legend.title          = element_text(size = 16, face = "bold"),
      legend.text           = element_text(size = 14),
      legend.box.spacing    = unit(0, "pt"),
      legend.margin         = margin(0, 0, 0, 0),
      legend.box.margin     = margin(-10, 0, 0, 0),
      legend.key.height     = unit(10, "pt"),
      legend.key.width      = unit(14, "pt"),
      
      # --- Plot margins (tight) ---
      plot.margin           = margin(0, 0, 0, 0)) + 
    guides(color = guide_legend(nrow = 2))
  
  ggsave(file = paste0('./output/plots/tern_plot_', vars[i], '_habitat.png'), 
         width = 5, 
         height = 5,
         dpi = 600)
  
  if(vars2[i] %in% c('Ammonia', 'Phosphate', 'Nitrite plus Nitrate')){
    readr::write_rds(m, paste0('./data/combined_figures/tern_', vars2[i], '.rds'))
  }
}

