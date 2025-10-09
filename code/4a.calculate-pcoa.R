library(tidyverse)
library(phyloseq)
library(ggplot2)
library(RColorBrewer)

level <- 'asv'

physeq <- 
  readr::read_rds(paste0('./data/ati-physeq-nonrrfy-', level, '.rds'))

#Make into relative abundance to run ordinations
physeq.ra <- transform_sample_counts(physeq, function(x) x/ sum(x))
data.ra <- as(sample_data(physeq.ra), "data.frame")

##Run a bray curtis pcoa and plot
bc.pcoa <- phyloseq::ordinate(physeq.ra, "PCoA", "bray")
readr::write_rds(bc.pcoa, 
                 paste0('./data/ati-pcoa-bray-object-', level, '.rds'))

df <- data.frame(data.ra,
                 PCoA1 = bc.pcoa$vectors[,1],
                 PCoA2 = bc.pcoa$vectors[,2],
                 PCoA3 = bc.pcoa$vectors[,3])

readr::write_rds(df, 
                 paste0('./data/ati-pcoa-bray-3-axes-metadata-', level, '.rds'))

#Run a unifrac curtis pcoa and plot
uni.pcoa <- phyloseq::ordinate(physeq.ra, "PCoA", "unifrac", weighted = TRUE)
readr::write_rds(uni.pcoa, 
                 paste0('./data/ati-pcoa-unifrac-object-', level, '.rds'))

df <- data.frame(data.ra,
                 PCoA1 = uni.pcoa$vectors[,1],
                 PCoA2 = uni.pcoa$vectors[,2],
                 PCoA3 = uni.pcoa$vectors[,3])

readr::write_rds(df, 
                 paste0('./data/ati-pcoa-unifrac-3-axes-metadata-', level, '.rds'))





