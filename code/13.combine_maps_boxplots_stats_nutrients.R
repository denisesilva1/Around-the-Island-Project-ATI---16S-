library(ggplot2)
library(patchwork)


# Ammonia -----------------------------------------------------------------
ammonia.maps <- 
  readr::read_rds('./data/combined_figures/Ammonia.rds') #+

ammonia.boxplots <-
  readr::read_rds('./data/combined_figures/ammonia_boxplots.rds')

ammonia <- ammonia.maps /
  ((plot_spacer() | ammonia.boxplots | plot_spacer()) +
     plot_layout(widths = c(.1, 10, .55)))

ammonia

ggsave(plot = ammonia, 
       filename = './output/maps/ammonia_maps_bxplt_21.22.23.png',
       dpi = 1000,
       width = 11,
       height = 7)

ggsave(plot = ammonia, 
       filename = './output/maps/ammonia_maps_bxplt_21.22.23.pdf',
       width = 11,
       height = 7)

# Nitrate+Nitrite ---------------------------------------------------------
nn.maps <- 
  readr::read_rds('./data/combined_figures/Nitrite_plus_Nitrate.rds') #+

nn.boxplots <-
  readr::read_rds('./data/combined_figures/nitrate_nitrite_boxplots.rds')

nn <- nn.maps /
  ((plot_spacer() | nn.boxplots | plot_spacer()) +
     plot_layout(widths = c(.1, 10, .55)))

nn

ggsave(plot = nn, 
       filename = './output/maps/nitrate_nitrite_maps_bxplt_21.22.23.png',
       dpi = 1000,
       width = 11,
       height = 7)

ggsave(plot = nn, 
       filename = './output/maps/nitrate_nitrite_maps_bxplt_21.22.23.pdf',
       width = 11,
       height = 7)

# Phosphate ---------------------------------------------------------------

phosph.maps <- 
  readr::read_rds('./data/combined_figures/Phosphate.rds') #+

phosph.boxplots <-
  readr::read_rds('./data/combined_figures/phosphate_boxplots.rds')

phosph <- phosph.maps /
  ((plot_spacer() | phosph.boxplots | plot_spacer()) +
     plot_layout(widths = c(.1, 10, .55)))

phosph

ggsave(plot = phosph, 
       filename = './output/maps/phosphate_maps_bxplt_21.22.23.png',
       dpi = 1000,
       width = 11,
       height = 7)

ggsave(plot = phosph, 
       filename = './output/maps/phosphate_maps_bxplt_21.22.23.pdf',
       width = 11,
       height = 7)


