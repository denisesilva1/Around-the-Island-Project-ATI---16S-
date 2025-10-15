library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(tidyr)
library(dplyr)
library(colorspace)

topn <- 30
level <- 'family'; level2 <- 'Family'

pal <- brewer.pal(12,"Paired") 
pal <- colorRampPalette(pal)(topn)

physeq_nr <- readr::read_rds(paste0('./data/ati-physeq-nonrrfy-', level, '.rds')) 

# Transform to relative abundance
phyloseq_rel_abund <- transform_sample_counts(physeq_nr, function(x) (x / sum(x))*100)

# Identify the top N taxa
top_taxa <- names(sort(taxa_sums(phyloseq_rel_abund), decreasing = TRUE)[1:topn])
readr::write_rds(top_taxa, paste0('./data/top-', topn,'-', level, '.rds'))

# Prune to keep only the top N genera
phyloseq_top <- prune_taxa(top_taxa, phyloseq_rel_abund)

# Prepare data for plotting
plot_data <- psmelt(phyloseq_top) %>%
  group_by(Habitat, Year, !!sym(level2)) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()

# Create the barplot
ggplot(plot_data, aes(x = Habitat, y = Abundance, fill = !!sym(level2))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pal) +
  labs(title = paste0("Relative Abundance of Top ", topn, " ", level2, " by Habitat"),
       x = "Habitat",
       y = "Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggh4x::facet_wrap2(~ Year, ncol = 1,
                     axes = "y", remove_labels = "y",
                     scale = "free_y") +
  coord_flip() +
  ylab("Mean Relative Abundance (%)") + 
  xlab("Habitat") + 
  theme_bw() + 
  theme(legend.position = 'bottom',
        axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold"))

ggplot2::ggsave(paste0('./output/plots/ati_relabund_', level, '_top', topn, '_21.22.23.png'), 
                height = 400, 
                width = 850,
                units = "mm",
                scale = 0.5, 
                dpi = 1000)

