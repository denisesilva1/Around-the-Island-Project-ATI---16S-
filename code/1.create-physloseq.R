library(tidyverse)
library(phyloseq)
library(decontam)
library(microbiome)
library(ggplot2)
library(qiime2R)
library(microViz)

## Aggregate/write metadata --------------------------------------------------

### Mismatched IDs
correct <- readr::read_csv('./data/corrected-sample-ids.csv')

### 2023 metadata
meta.columns <- 
  c('sample-id', 'Sample_type', 'Site', 'UniqueID', 'Year', 'TimePeriod', 
    'Phosphate', 'Silicate', 'Nitrite_plus_Nitrate', 'Ammonia', 'Percent_C', 
    'Percent_N', 'C_to_N_ratio', 'Latitude', 'Longitude', 'Island_shore', 'Habitat', 
    'Lagoon', 'Distance_to_crest', 'Distance_to_shore', 'Distance_to_pass',
    'Distance_to_deep_lagoon_water', 'Distance_to_population_center', 'Population_center')

turb <- 
  readr::read_csv('./data/metadata/Turb_CHN_compiled_April2023.csv') %>% 
  dplyr::select(-Sample_ID, -Weight_ug, -Percent_H) %>% 
  mutate(Site_Number = case_when(
    Site_Number == 'Cooks 1' ~ 'Cook 1',
  #   Site_Number=='OPU 1' ~ 'Opu1',
  #   Site_Number=='OPU 2' ~ 'Opu2',
    Site_Number=='LTER 3 Back Reef' ~ 'LTER_3_Backreef',
    Site_Number=='LTER 3 Fringe Reef' ~ 'LTER_3_Fringe',
    Site_Number=='LTER 6 Fringe' ~ 'LTER_6_Fringe',
    TRUE ~ Site_Number
  )) %>%
  rename(Site = Site_Number)

dN15 <- readr::read_csv('./data/metadata/ati-21.22.23-dn15-dc13.csv') %>% 
  #dplyr::filter(Year == 2023) %>% 
  dplyr::select(Site, dN15, dC13) %>% 
  na.omit() %>% 
  dplyr::filter(!duplicated(Site))

ati_meta <- list()

ati_meta[[1]] <- 
  readr::read_delim('./data/metadata/metadata_ATI_2023.txt') %>% 
  dplyr::select(-Percent_C, -Percent_N, -C_to_N_ratio) %>% 
  dplyr::rename(Site = Sites,
                Nitrite_plus_Nitrate = `Nitrite+Nitrate`) %>% 
  dplyr::filter(!Site == 'Positive control') %>% 
  dplyr::mutate(Latitude = ifelse(Site == 'GUMP', -17.4909, Latitude),
                Longitude = ifelse(Site == 'GUMP', -149.826, Longitude),
                Sample_type = ifelse(Site == 'Negative control', 
                                     'control',
                                     'lagoon water'),
                Site = case_when(
                  Site == 'Cooks1' ~ 'Cook 1',
                  Site== 'Opu1' ~ 'OPU 1',
                  Site=='Opu2' ~ 'OPU 2',
                #   Site=='LTER 3 Back Reef' ~ ,
                #   Site=='LTER_3_Backreef'~'LTER 3 Fringe Reef',
                #   Site=='LTER_6_FRINGE'~'LTER 6 Fringe',
                  TRUE ~ Site)
                ) %>% 
  dplyr::left_join(turb, join_by(Site)) %>% 
  dplyr::rename('sample-id' = '#SampleID' ) %>% 
  dplyr::mutate(UniqueID = paste0(`sample-id`, '_', TimePeriod),
                `sample-id` = paste0(`sample-id`, '_23')) %>% 
  dplyr::select(all_of(meta.columns)) %>% 
  dplyr::left_join(dN15, join_by(Site))

### 2021/2022 metadata
years <- c(2021, 2022)
yrs <- c('_21', '_22')

meta_https <- 
  'https://raw.githubusercontent.com/njsilbiger/ATI_NutrientRegimes/main/Data/AllNutrientData_clusters.csv'

meta_files <- list.files(path = './data/metadata',
                         pattern = '2021|2022', 
                         full.names = TRUE) 

for(i in 1:length(meta_files)){
  
  dN15 <- readr::read_csv('./data/metadata/ati-21.22.23-dn15-dc13.csv') %>% 
    dplyr::filter(Year == years[i]) %>% 
    dplyr::select(Site, dN15, dC13) %>% 
    na.omit() %>% 
    dplyr::filter(!duplicated(Site))
  
  ati_meta[[i+1]] <- 
    readr::read_csv(meta_https) %>% 
    dplyr::filter(Year == years[i]) %>% 
    dplyr::full_join(
      readr::read_delim(meta_files[i]),
      dplyr::join_by(Site)
    ) %>% 
    #dplyr::rename('sample-id'  = 'Sample-id')
    {
      if (!'sample-id' %in% colnames(.)) {
        dplyr::rename(., 'sample-id' = 'Sample-id')
      } else {
        .
      }
    } %>% 
    mutate(`sample-id` = ifelse(is.na(`sample-id`), Site, `sample-id`)) %>% 
    relocate('sample-id') %>% 
    dplyr::mutate(UniqueID = paste0(`sample-id`, '_', TimePeriod),
                  `sample-id` = paste0(`sample-id`, yrs[i])) %>% 
    dplyr::select(all_of(meta.columns)) %>% 
    dplyr::left_join(dN15, join_by(Site))
}

ati_meta <- 
  ati_meta %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(Sample_type = ifelse(Sample_type == 'lagoon_water',
                                     'lagoon water', Sample_type)) %>% 
  dplyr::filter(!Site %in% c('CG1', 'GUMP', 'LTER_4_Backreef'))

# Replace IDs using match
ati_meta$`sample-id` <- 
  ifelse(ati_meta$`sample-id` %in% correct$meta_id, 
         correct$qza_id[match(ati_meta$`sample-id`, correct$meta_id)], 
         ati_meta$`sample-id`)

ati_meta$Sample_type[ati_meta$`sample-id` == '56_22'] <- 'lagoon water'

# tmp <- ati_meta %>% 
#   dplyr::select(`sample-id`, Site, UniqueID, Latitude, Longitude)
# 
# write.csv(tmp, 'C:/workfolder/ati-lat-lon.csv', row.names = F)

# ati_meta %>% 
#   dplyr::filter(`sample-id` %in% correct$qza_id) %>% 
#   View()

# ati_meta <- 
#   ati_meta %>% 
#   dplyr::select(`sample-id`, Site, UniqueID, Habitat, Sample_type)

readr::write_tsv(ati_meta, './data/metadata/ati-21.22.23-metadata.tsv')

## Generate phyloseq object -------------------------------------------------
features <- './data/qiime2-pipeline/ATI_21_22_23-filtered_table.qza'
tree <- './data/qiime2-pipeline/ATI_21_22_23-rooted-tree.qza'
taxonomy <- './data/qiime2-pipeline/ATI_21_22_23-tax.qza'
metadata <- './data/metadata/ati-21.22.23-metadata.tsv'

# test <- qiime2R::qza_to_phyloseq(features = features) %>% 
#   colnames() %>% 
#   data.frame() 
# write.csv(test, 'qza_names.csv', row.names = F)
# 
# tmp <- ati_meta %>% dplyr::select(`sample-id`)
# write.csv(tmp, 'meta_names.csv', row.names = F)

physeq <- 
  qiime2R::qza_to_phyloseq(features = features, 
                           tree = tree, 
                           taxonomy = taxonomy, 
                           metadata = metadata) 

microViz::phyloseq_validate(physeq)

unknown_labels <- c("", "NA", "(?i).*uncultured.*", "Incertae Sedis",
                    "(?i).*Incertae.*", "uncultured bacterium",
                    "uncultured soil bacterium")

physeq <- 
  tax_fix(physeq, 
          unknowns = unknown_labels)

microViz::phyloseq_validate(physeq)

#4ef10f3fb5c4348a4e3a32a991a08904


# dat <- 
#   sample_data(physeq) %>% 
#   data.frame() %>%
#   mutate(`sample-id` = row.names(.)) %>% 
#   dplyr::select(`sample-id`, Site, UniqueID)
# 
# tmp <- ati_meta %>% 
#   full_join(dat,
#             join_by(UniqueID)) %>% 
#   dplyr::filter(is.na(Site.y))
# 
# write.csv(tmp, './data/mismatch-sample-ids.csv', row.names = FALSE)

### Remove unassigned reads (unassigned at the kingdom level), euks, chloroplast
physeq <- 
  physeq %>% 
  phyloseq::subset_taxa(!(Kingdom == 'Unassigned' | 
                            Kingdom == 'Eukaryota' | 
                            Order == 'Chloroplast' |
                            Family == 'Mitochondria' |
                            Species == 'Bacteria Kingdom'))



### Check the taxonomic classification is correct
phyloseq::rank_names(physeq) #7 ranks

## Decontamination ---------------------------------------------------------
### make a sample data frame
sample.data <- 
  as(phyloseq::sample_data(physeq), 'data.frame') %>% 
  dplyr::mutate(library_size = sample_sums(physeq)) %>% 
  dplyr::arrange(library_size) %>% 
  dplyr::mutate(Index = seq(nrow(.)),
                Habitat_type = ifelse(is.na(Habitat),
                                      Sample_type,
                                      Habitat)) 

### visualize library size
ggplot(sample.data) +
  geom_point(aes(x=Index, 
                 y=library_size, 
                 color = Sample_type),
             #alpha = 0.75,
             size=3) +
  scale_color_manual(values = c("control" = "black", 
                                "lagoon water" = "lightblue")) +
  theme_bw()

ggplot(sample.data %>% 
         filter(Index < 50)) +
  geom_point(aes(x=Index, 
                 y=library_size, 
                 color = Habitat_type),
             #alpha = 0.5,
             size=3) +
  scale_color_brewer(palette = 'Set1') +
  theme_bw()

### Check for contaminants
phyloseq::sample_data(physeq)$is.neg <- 
  phyloseq::sample_data(physeq)$Sample_type == 'control'

### prevalence and threshold 0.5 (more conservative)
contamdf.prev <- 
  decontam::isContaminant(physeq, 
                          method = 'prevalence', 
                          neg = 'is.neg', 
                          threshold = 0.5)

table(contamdf.prev$contaminant) # Number of TRUE potential contaminants

### Remove contaminants and negative controls
keep_taxa <- !contamdf.prev$contaminant[phyloseq::taxa_names(physeq)]
keep_taxa[is.na(keep_taxa)] <- TRUE
physeq <- phyloseq::prune_taxa(keep_taxa, physeq) %>% 
  phyloseq::subset_samples(!is.neg)

# physeq <- 
#   phyloseq::prune_taxa(!contamdf.prev$contaminant, physeq) %>% 
#   phyloseq::subset_samples(!is.neg)

## Additional cleaning and save -------------------------------------------

### Remove singletons
physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq) > 1, 
                               physeq)

#microbiome::summarize_phyloseq(physeq)

readr::write_rds(physeq, 
                 file = './data/ati-physeq-nonrrfy-asv.rds', 
                 compress = 'xz')

# physeq_nr <- physeq %>% 
#   phyloseq::tax_glom(taxrank = 'Species')
# 
# readr::write_rds(physeq_nr, 
#                  file = './data/ati-physeq-nonrrfy-species.rds', 
#                  compress = 'xz')
# 
# physeq_nr <- physeq %>% 
#   phyloseq::tax_glom(taxrank = 'Genus')
# 
# readr::write_rds(physeq_nr, 
#                  file = './data/ati-physeq-nonrrfy-genus.rds', 
#                  compress = 'xz')

physeq_nr <- physeq %>% 
  phyloseq::tax_glom(taxrank = 'Family')

readr::write_rds(physeq_nr, 
                 file = './data/ati-physeq-nonrrfy-family.rds', 
                 compress = 'xz')

### Make a rarefied phyloseq object for alpha diversity analyses
physeq <- 
  phyloseq::rarefy_even_depth(physeq, 
                              #sample.size = 10500, 
                              sample.size = 5000, 
                              rngseed = 102086) 

readr::write_rds(physeq, 
                 file = './data/ati-physeq-rrfy-asv.rds', 
                 compress = 'xz')

# physeq_rr <- physeq %>% 
#   phyloseq::tax_glom(taxrank = 'Species')
# 
# write_rds(physeq_rr, './data/ati-physeq-rrfy-species.rds', 
#           compress = 'xz')
# 
# physeq_rr <- physeq %>% 
#   phyloseq::tax_glom(taxrank = 'Genus')
# 
# write_rds(physeq_rr, './data/ati-physeq-rrfy-genus.rds', 
#           compress = 'xz')

physeq_rr <- physeq %>% 
  phyloseq::tax_glom(taxrank = 'Family')

write_rds(physeq_rr, './data/ati-physeq-rrfy-family.rds', 
          compress = 'xz')

