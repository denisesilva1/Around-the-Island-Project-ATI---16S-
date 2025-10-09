library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(tidyr)
library(dplyr)
library(spmodel)

level <- 'family'; level2 <- 'Family'

physeq_rrfy <- 
  readr::read_rds(paste0('./data/ati-physeq-rrfy-', level2, '.rds')) %>% 
  #transform_sample_counts(function(x) (x / sum(x))*100) %>% 
  subset_taxa(Family == "Vibrionaceae" | 
                Family == "Nitrincolaceae" | 
                Family == "Pseudoalteromonadaceae" | 
                Family == "Moraxellaceae" |
                Family == "Halomonadaceae")

taxa <- tax_table(physeq_rrfy)@.Data %>% 
  data.frame() %>% 
  dplyr::mutate(asv = row.names(.)) %>% 
  dplyr::select(asv, Family)

abundances <- 
  phyloseq::sample_data(physeq_rrfy) %>% 
  as_tibble() %>% 
  dplyr::select(Site, Latitude, Longitude, Year, Habitat, Island_shore,
                Silicate, Percent_N, Phosphate, Ammonia, Nitrite_plus_Nitrate,
                dN15, dC13) %>%
  dplyr::mutate(Site = as.factor(Site),
                Year = as.factor(Year),
                Habitat = factor(Habitat,
                                    levels = c('Fringing reef','Bay','Mid lagoon','Reef crest','Reef pass'))) %>% 
  bind_cols(
    phyloseq::otu_table(physeq_rrfy) %>% 
      t()) %>% 
  na.omit() %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(crs = 3305)

# Verify taxa names and order
#colnames(abundances)[4:ncol(abundances)] 
colnames(abundances)[grep(paste(taxa$asv, collapse = "|"),
                          colnames(abundances))] <- taxa$Family


# Vibrionaceae ---------------------------------------------------------------------
formula <- 
  log10(Vibrionaceae + 1) ~
  Habitat +
  Year +
  Silicate + Percent_N + Phosphate + Ammonia + Nitrite_plus_Nitrate

sp.mod0 <- splm(formula = formula, 
                data = abundances,
                #partition_factor = ~ Island_shore,
                random = ~ Site,
                #estmethod='ml',
                spcov_type = 'gaussian') 

#sp.mod1 <- update(sp.mod0, . ~ . + Phosphate )

glances(sp.mod0)
summary(sp.mod0)
loocv(sp.mod0)

# Nitrincolaceae ---------------------------------------------------------------------
formula <- update(formula, log10(Nitrincolaceae + 1) ~ .)

sp.mod0 <- splm(formula = formula, 
                data = abundances,
                #partition_factor = ~ Island_shore,
                random = ~ Site,
                #estmethod='ml',
                spcov_type = 'gaussian') 

glances(sp.mod0)
summary(sp.mod0)
loocv(sp.mod0)

# Pseudoalteromonadaceae ---------------------------------------------------------------------
formula <- update(formula, log10(Pseudoalteromonadaceae + 1) ~ .)

sp.mod0 <- splm(formula = formula, 
                data = abundances,
                #partition_factor = ~ Island_shore,
                random = ~ Site,
                #estmethod='ml',
                spcov_type = 'gaussian') 

glances(sp.mod0)
summary(sp.mod0)
loocv(sp.mod0)

# Moraxellaceae ---------------------------------------------------------------------
formula <- update(formula, log10(Moraxellaceae + 1) ~ .)

sp.mod0 <- splm(formula = formula, 
                data = abundances,
                #partition_factor = ~ Island_shore,
                random = ~ Site,
                #estmethod='ml',
                spcov_type = 'gaussian') 

glances(sp.mod0)
summary(sp.mod0)
loocv(sp.mod0)

# Halomonadaceae ---------------------------------------------------------------------
formula <- update(formula, log10(Halomonadaceae + 1) ~ .)

sp.mod0 <- splm(formula = formula, 
                data = abundances,
                #partition_factor = ~ Island_shore,
                random = ~ Site,
                #estmethod='ml',
                spcov_type = 'gaussian') 

glances(sp.mod0)
summary(sp.mod0)
loocv(sp.mod0)

