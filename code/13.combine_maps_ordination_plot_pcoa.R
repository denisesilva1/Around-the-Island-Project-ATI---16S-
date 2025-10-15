library(dplyr)
library(ggplot2)
library(patchwork)

pcoa <- readr::read_rds('./data/combined_figures/combined_pcoa_habitat_by_year.rds')

maps <- readr::read_rds('./data/combined_figures/combined_pcoa_rgb_map_by_year.rds')

combined <- 
  maps / 
  ((plot_spacer() | pcoa | plot_spacer()) +
     plot_layout(widths = c(-.3, 20, .87)))

ggsave(plot = combined, 
       filename = './output/maps/pcoa_map_plots_21.22.23b.png',
       dpi = 1000,
       width = 11,
       height = 7)
