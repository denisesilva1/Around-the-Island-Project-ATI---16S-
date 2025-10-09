library(phyloseq)
library(ggplot2)
library(tidyverse)
library(picante)
library(vegan)


# Data prep ---------------------------------------------------------------
level <- 'asv'

physeq_rrfy <-
  readr::read_rds(paste0('./data/ati-physeq-rrfy-', level, '.rds'))

data.rrfy <- as(sample_data(physeq_rrfy), "data.frame") 

OTU <- 
  as(phyloseq::otu_table(physeq_rrfy), 
     "matrix") %>% 
  t()

# Calculate alpha diversity metrics ---------------------------------------
rich <- 
  vegan::estimateR(OTU) %>% 
  t() %>% 
  data.frame() %>% 
  dplyr::select(S.obs, S.chao1) 

faith.pd <- 
  picante::pd(samp = OTU,
              tree = phyloseq::phy_tree(physeq_rrfy), 
              include.root = FALSE) %>% 
  dplyr::select(PD)

shannon <- 
  vegan::diversity(OTU, 
            index = "shannon") %>% 
  as_tibble() %>% 
  rename(shannon = value)


# Combine alpha diversity metrics -----------------------------------------
lookup <- c('Microbial_Richness'='S.obs',
            'Microbial_Chao1' = 'S.chao1',
            'Microbial_Shannon_Diversity'='shannon',
            'Microbial_Phylogenetic_Diversity'='PD',
            'Microbial_Evenness'='evenness')

alpha_div <- 
  cbind(rich, shannon, faith.pd) %>% 
  dplyr::mutate(evenness = shannon / log(S.obs)) %>% 
  cbind(data.rrfy) %>% 
  dplyr::rename(all_of(lookup)) %>% 
  dplyr::relocate(Site, Year) %>% 
  dplyr::filter(!is.na(Latitude)) %>% 
  as_tibble()

write_rds(alpha_div, 
          file = paste0('./data/ati-alpha-diversity-metrics-', level, '.rds'), 
          compress = "xz")
