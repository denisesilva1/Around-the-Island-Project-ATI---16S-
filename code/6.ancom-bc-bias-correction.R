library(phyloseq)
library(tidyverse)
library(ANCOMBC)
library(microbiome)
library(lme4)

# Data prep ---------------------------------------------------------------
topn <- 20
#level <- 'asv'; level2 <- 'ASV'
#level <- 'species'; level2 <- 'Species'
#level <- 'genus'; level2 <- 'Genus'
level <- 'family'; level2 <- 'Family'

physeq_nr <- 
  readr::read_rds(paste0('./data/ati-physeq-nonrrfy-', level, '.rds'))

# taxa <- tax_table(physeq_nr)
# View(taxa)

# Check for correlated taxa
otu_table_raw <-
  as(otu_table(physeq_nr), 'matrix')
cor_matrix <-
  cor(otu_table_raw, method = 'pearson') # or 'spearman'

ggplot(data = cor_matrix[upper.tri(cor_matrix)] %>%
         data.frame()) +
  geom_histogram(aes(x = .)) +
  ggtitle('Histogram of Taxa Correlations') +
  xlab('Correlations') +
  theme_bw()

metadata <- 
  sample_data(physeq_nr) %>% 
  data.frame() %>% #names()
  dplyr::select(Site, Year, Habitat, Island_shore) %>% 
  dplyr::mutate(samp_id = row.names(.),
                Year = factor(Year),
                Habitat = factor(Habitat),
                Island_shore = factor(Island_shore)) %>% 
  dplyr::filter(!is.na(Habitat))

metadata$Site[metadata$Site == 'Opu1'] <- c('OPU 1', 'OPU 1')
metadata$Site[metadata$Site == 'Opu2'] <- 'OPU 2'
metadata$Site[metadata$Site == 'Cooks1'] <- 'Cook 1'

# Count observations per SiteID
site_counts <- table(metadata$Site)

valid_sites <-
  names(site_counts[site_counts > 1])
metadata <-
  metadata[metadata$Site %in% valid_sites, ]

# Subset the phyloseq object
physeq_nr <-
  prune_samples(sample_names(physeq_nr) %in% rownames(metadata), physeq_nr)

sample_data(physeq_nr) <- metadata

habitats <- metadata %>% 
  pull(Habitat) %>%
  unique() %>% 
  as.character()

# Define new control settings
control <- lmerControl(
  optimizer = "nloptwrap",
  optCtrl = list(method = "nlminb",
                 maxfun = 1e6))

control <- 
  lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))

for(i in 1:length(habitats)){
  print(habitats[i])
  
  phsq_tmp <- physeq_nr %>%
    subset_samples(Habitat == habitats[i])
  
  ancombc2_results <-
    ancombc2(data = phsq_tmp,
             tax_level = level2,
             fix_formula = 'Year',
             #rand_formula = '(1 | Site)',
             #lme_control = control,
             p_adj_method = 'BY', #Many correlated taxa
             group = 'Year',
             struc_zero = FALSE,
             neg_lb = TRUE,
             alpha = 0.05,
             prv_cut = 0.1,
             n_cl = 8,
             verbose = TRUE,
             pairwise = FALSE)
  
  readr::write_rds(ancombc2_results,
                   paste0('./data/ancombc/ancombc2-',
                          habitats[i],'-pairwise-Year-', level, '.rds'),
                   compress = 'xz')
  
}

ancombc_results <- 
  list.files(path = './data/ancombc/', 
             pattern = level,
             full.names = TRUE) %>% 
  map_dfr(~ {
    # Read the .rds file
    file_data <- readr::read_rds(.x)
    # Extract the "res" component and return it
    res_data <- file_data$res
    # Add a column with the file name (or a custom identifier)
    habitat_name <- .x %>%
      basename() %>%  # Extract the file name without the path
      str_remove("^ancombc2-") %>%    # Remove the prefix "ancombc2-"
      str_remove(paste0("-pairwise-Year-", level, "\\.rds$")) # Remove the suffix "-pairwise-Year-genus.rds"
    # Add the habitat name as a column
    res_data %>%
      mutate(Habitat = habitat_name)
  })


# Reshape the results into long format
results_long <- 
  ancombc_results %>%
  dplyr::select(taxon, lfc_Year2022, 
                lfc_Year2023, q_Year2022, 
                q_Year2023, Habitat) %>% 
  pivot_longer(cols = starts_with("lfc_") | starts_with("q_"),
               names_to = c(".value", "Year"), 
               names_sep = "_") %>%
  dplyr::mutate(significant = q < 0.05,
                Year = str_remove(Year, patter = "Year")) 

if(level == 'genus'){
  toptax <- 21
}else if(level == 'species'){
  toptax <- 26
}else{
  toptax <- 60
}

lfc <- results_long %>% 
  group_by(taxon) %>% 
  summarise(effect_size = mean(abs(lfc))) %>% 
  dplyr::arrange(desc(abs(effect_size))) %>% 
  dplyr::mutate(top_effect_size = case_when(
    row_number() <= toptax ~ TRUE,
    TRUE ~ FALSE
  ))

results_long <- results_long %>% 
  dplyr::left_join(lfc,
                   join_by(taxon)) %>% 
  #dplyr::arrange(desc(significant), desc(abs(effect_size))) %>% 
  dplyr::mutate(taxon2 = case_when(
    top_effect_size == TRUE & significant == TRUE ~ taxon, 
    significant == FALSE ~ "Non-Significant", 
    TRUE ~ "Significant, Other" 
  ))

if(level == 'genus'){
  results_long$taxon2[grep('Alphaproteobacteria_Incertae Sedis', results_long$taxon2)] <- 'Alphaproteobacteria'
  results_long$taxon2[grep('Magnetospiraceae', results_long$taxon2)] <- 'Magnetospiraceae'
  #results_long$taxon2[grep('Defluviicoccales', results_long$taxon2)] <- 'Unclassified, Order - Defluviicoccales'
  results_long$taxon2[grep('PS1 clade_Incertae Sedis', results_long$taxon2)] <- 'PS1 clade'
  #results_long$taxon2[grep('Pirellulaceae_Incertae', results_long$taxon2)] <- 'Unclassified, Family - Pirellulaceae'
  results_long$taxon2[grep('Verrucomicrobiales_DEV007', results_long$taxon2)] <- 'DEV007'
  results_long$taxon2[grep('Dehalococcoidia_SAR202', results_long$taxon2)] <- 'SAR202 clade'
  results_long$taxon2[grep('Parvibaculales_PS1', results_long$taxon2)] <- 'PS1 clade'
  results_long$taxon2[grep('Flavobacteriaceae_Incertae Sedis', results_long$taxon2)] <- 'Flavobacteriaceae'
}

if(level == 'species'){
  results_long$taxon2[grep('Marinoscillum_uncultured bacterium', results_long$taxon2)] <- 'Marinoscillum'
  results_long$taxon2[grep('NS4 marine group_uncultured Flavobacteriaceae', results_long$taxon2)] <- 'Flavobacteriaceae'
  results_long$taxon2[grep('Rubripirellula_uncultured bacterium', results_long$taxon2)] <- 'Rubripirellula'
  results_long$taxon2[grep('Alphaproteobacteria_Incertae Sedis_Incertae', results_long$taxon2)] <- 'Alphaproteobacteria'
  results_long$taxon2[grep('Alphaproteobacteria_Parvibaculales_PS1 clade_Incertae', results_long$taxon2)] <- 'Parvibaculales PS1 clade'
  results_long$taxon2[grep('Candidatus Endoecteinascidia_uncultured gamma proteobacterium', results_long$taxon2)] <- 'Candidatus'
  results_long$taxon2[grep('DEV007_Incertae Sedis', results_long$taxon2)] <- 'DEV007'
  results_long$taxon2[grep('uncultured Sphingobacteriia bacterium', results_long$taxon2)] <- 'Sphingobacteriia'
  results_long$taxon2[grep('Flavobacteriaceae_NS4', results_long$taxon2)] <- 'Flavobacteriia marine group'
  results_long$taxon2[grep('uncultured Chloroflexi bacterium', results_long$taxon2)] <- 'Chloroflexi'
  results_long$taxon2[grep('Phycisphaeraceae_Urania-1B-19', results_long$taxon2)] <- 'Urania-1B-19 marine sediment group'
  results_long$taxon2[grep('Pirellulales_Pirellulaceae', results_long$taxon2)] <- 'Pirellulaceae'
  results_long$taxon2[grep('Lentimonas_uncultured', results_long$taxon2)] <- ' Lentimonas'
  results_long$taxon2[grep('SCGC AAA164-E04_uncultured', results_long$taxon2)] <- 'SCGC AAA164-E04'
  results_long$taxon2[grep('SAR116 clade_Incertae', results_long$taxon2)] <- 'SAR116 clade'
  results_long$taxon2[grep('Magnetospiraceae_Incertae', results_long$taxon2)] <- 'Magnetospiraceae'
}

if(level == 'family'){
  results_long$taxon2[grep('Bacteroidia_Chitinophagales', results_long$taxon2)] <- 'Chitinophagales'
  results_long$taxon2[grep('Dehalococcoidia_SAR202 clade', results_long$taxon2)] <- 'SAR202 clade'
  results_long$taxon2[grep('Alphaproteobacteria_Incertae', results_long$taxon2)] <- 'SAR202 clade'
  results_long$taxon2[grep('Planctomycetota_vadinHA49', results_long$taxon2)] <- 'vadinHA49'
  results_long$taxon2[grep('Alphaproteobacteria_Defluviicoccales', results_long$taxon2)] <- 'Defluviicoccales'
  
}

tax <- results_long %>% 
  dplyr::filter(!taxon2 %in% c('Significant, Other', 'Non-Significant')) %>% 
  pull(taxon2) %>% 
  unique()

results_long <- results_long %>% 
  mutate(taxon2 = factor(taxon2,
                         levels = c(tax, c('Significant, Other', 'Non-Significant'))))

library(RColorBrewer)
pal <- brewer.pal(12,"Paired")
pal <- colorRampPalette(pal)(length(unique(results_long$taxon2))-2)
pal <- colorspace::darken(pal, amount = 0.2)
pal <- c(pal, c("#b5eaff", "#bdbdbd"))
# 
# pal <- brewer.pal(9, "Set1") %>% # 9 colors from Set1
#   c(brewer.pal(8, "Set2"), brewer.pal(1, "Set3")[1]) 
# pal <- colorspace::darken(pal, amount = 0.3)
#   
# pal <- c(pal, c("#b5eaff", "lightgrey"))

# Volcano Plots with Facets
ggplot(results_long %>% 
         dplyr::filter(!taxon2 %in% c('Significant, Other', 'Non-Significant')), 
       aes(x = lfc, 
           y = -log10(q), 
           color = taxon2)) +
  geom_point(size = 3, alpha = 0.9, stroke = NA) +
  facet_grid(rows = vars(Habitat), cols = vars(Year)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Significance threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +      # Log-fold change thresholds
  # facet_wrap(~ Year, 
  #            #scales = "free_y", 
  #            ncol = 2) +                         # Facet by variable
  #scale_color_manual(values = c("#bdbdbd", "#1d91c0")) +                               # Gray = not significant, Red = significant
  scale_color_manual(values = pal) +
  labs(
    #title = "Volcano Plots for Pairwise Comparisons Across Variables",
    x = "Log-Fold Change",
    y = "-log10 Adjusted p-value"
  ) +
  theme_bw() +
  ggpubr::labs_pubr() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    strip.background = element_blank(),
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(ncol = 4))

ggsave(paste0('./output/plots/ancombc2-volcano-year-', level, '-top35.png'),
       dpi = 600,
       units = 'mm', 
       width = 230,
       height = 250)




top <- results_long %>% 
  filter(variable %in% c(2022, 2023),
         lfc > 1 | lfc < -1,
         significant == TRUE) #%>% 
# group_by(taxon) %>% 
# summarise(abs_lfc = mean(abs(lfc))) %>% 
# #mutate(abs_lfc = abs(lfc)) %>%  # Compute absolute LFC
# arrange(desc(abs_lfc)) %>%      # Sort by absolute LFC in descending order
# slice_head(n = 10)              # Select the top 10 rows

results_year <- results_long %>% 
  filter(variable %in% c(2022, 2023)) %>% 
  mutate(taxon2 = ifelse(taxon %in% top$taxon, taxon, 'Other'))

if(level == 'genus'){
  results_year$taxon2[grep('Alphaproteobacteria_Incertae Sedis', results_year$taxon2)] <- 'Unclassified, Class - Alphaproteobacteria'
  results_year$taxon2[grep('Magnetospiraceae', results_year$taxon2)] <- 'Unclassified, Family - Magnetospiraceae'
  #results_year$taxon2[grep('Defluviicoccales', results_year$taxon2)] <- 'Unclassified, Order - Defluviicoccales'
  results_year$taxon2[grep('PS1 clade_Incertae Sedis', results_year$taxon2)] <- 'Unclassified, Family - PS1 clade'
  #results_year$taxon2[grep('Pirellulaceae_Incertae', results_year$taxon2)] <- 'Unclassified, Family - Pirellulaceae'
  results_year$taxon2[grep('Verrucomicrobiales_DEV007', results_year$taxon2)] <- 'Unclassified, Family - DEV007'
}

pal <- RColorBrewer::brewer.pal(12,"Paired") 
pal <- colorRampPalette(pal)(nrow(top))

ggplot() +
  geom_point(data = results_year %>% 
               filter(taxon2 == 'Other'),
             aes(x = lfc,
                 y = -log10(q)),
             color = 'grey', 
             size = 3, alpha = 0.5, stroke = NA) +
  geom_point(data = results_year %>% 
               filter(taxon2 != 'Other'), 
             aes(x = lfc, 
                 y = -log10(q), 
                 color = taxon2),
             size = 5, alpha = 1, stroke = NA) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Significance threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +      # Log-fold change thresholds
  facet_wrap(~ variable, 
             #scales = "free_y", 
             ncol = 2) +                         # Facet by variable
  scale_color_manual(values = pal) +                               # Gray = not significant, Red = significant
  labs(
    #title = "Volcano Plots for Pairwise Comparisons Across Variables",
    x = "Log-Fold Change",
    y = "-log10 Adjusted p-value"
  ) +
  theme_minimal() +
  ggpubr::labs_pubr() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave(paste0('./output/plots/ancombc2-year-volcano-by-taxon-', level, '.png'),
       dpi = 600,
       height = 5,
       width = 11)

top <- results_long %>% 
  filter(!variable %in% c(2021, 2023)) %>% 
  mutate(abs_lfc = abs(lfc)) %>%  # Compute absolute LFC
  arrange(desc(abs_lfc)) %>%      # Sort by absolute LFC in descending order
  slice_head(n = 10)              # Select the top 10 rows

results_year <- results_long %>% 
  filter(variable %in% c(2021, 2023)) %>% 
  mutate(taxon2 = ifelse(taxon %in% top$taxon, taxon, 'Other'))

if(level == 'genus'){
  results_year$taxon[results_year$taxon == 'Bacteria_Pseudomonadota_Alphaproteobacteria_Incertae Sedis_Incertae Sedis_Incertae Sedis'] <-'Unclassified, Class - Alphaproteobacteria'
  results_year$taxon[grep('Magnetospiraceae', results_year$taxon)] <- 'Unclassified, Family - Magnetospiraceae'
}

pal <- RColorBrewer::brewer.pal(12,"Paired") 
pal <- colorRampPalette(pal)(nrow(top))

ggplot() +
  geom_point(data = results_year %>% 
               filter(taxon2 == 'Other'),
             aes(x = lfc,
                 y = -log10(q)),
             color = 'grey', 
             size = 3, alpha = 0.3, stroke = NA) +
  geom_point(data = results_year %>% 
               filter(taxon2 != 'Other'), 
             aes(x = lfc, 
                 y = -log10(q), 
                 color = taxon2),
             size = 3, alpha = 0.9, stroke = NA) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Significance threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +      # Log-fold change thresholds
  facet_wrap(~ variable, 
             #scales = "free_y", 
             ncol = 2) +                         # Facet by variable
  scale_color_manual(values = pal) +                               # Gray = not significant, Red = significant
  labs(
    #title = "Volcano Plots for Pairwise Comparisons Across Variables",
    x = "Log-Fold Change",
    y = "-log10 Adjusted p-value"
  ) +
  theme_minimal() +
  ggpubr::labs_pubr() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave(paste0('./output/plots/ancombc2-year-volcano-by-taxon-', level, '.png'),
       dpi = 600,
       height = 4,
       width = 10.2)


res <- ancombc2_results %>% 
  purrr::pluck('res')

filtered_df <- results[apply(results[, grep('diff', names(results))], 1, any), ]

# Generate Volcano Plot
ggplot(results, 
       aes(x = LFC, 
           y = -log10(p_value), 
           color = significant)) +
  geom_point() +
  geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
  labs(title = 'Volcano Plot of Habitat Comparisons',
       x = 'Log-Fold Change (LFC)',
       y = '-log10(p-value)') +
  theme_minimal()

