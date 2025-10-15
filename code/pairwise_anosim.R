# --- deps ---
library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

pairwise_anosim_by_habitat <- function(
    physeq,
    dist_method    = "bray",
    max_perm       = 999,      # cap permutations
    min_group_n    = 3,        # min samples per group to run ANOSIM
    seed           = 20252708, # RNG seed
    avoid_enum_msg = TRUE      # avoid "complete enumeration" note
){
  stopifnot(inherits(physeq, "phyloseq"))
  
  samp <- as(sample_data(physeq), "data.frame")
  need <- c("Year","Island_shore","Habitat")
  miss <- setdiff(need, colnames(samp))
  if (length(miss)) stop("Missing required sample vars: ", paste(miss, collapse=", "))
  
  # factor-ize
  samp$Year         <- as.factor(samp$Year)
  samp$Island_shore <- as.factor(samp$Island_shore)
  samp$Habitat      <- as.factor(samp$Habitat)
  sample_data(physeq)$Year         <- samp$Year
  sample_data(physeq)$Island_shore <- samp$Island_shore
  sample_data(physeq)$Habitat      <- samp$Habitat
  
  habitats <- levels(samp$Habitat)
  if (is.null(habitats)) habitats <- unique(samp$Habitat)
  
  # helper to avoid Inf from choose()
  safe_choose <- function(n, k) {
    val <- suppressWarnings(choose(n, k))
    if (is.infinite(val)) exp(lchoose(n, k)) else val
  }
  
  out_list <- list()
  k <- 1
  
  for (h in habitats) {
    # subset to this habitat WITHOUT NSE (no "h not found")
    keep_idx <- which(as.character(sample_data(physeq)$Habitat) == as.character(h))
    if (length(keep_idx) < 2) next
    phy_h <- prune_samples(sample_names(physeq)[keep_idx], physeq)
    phy_h <- prune_samples(sample_sums(phy_h) > 0, phy_h)
    if (nsamples(phy_h) < 2) next
    
    # distance
    dist_h <- phyloseq::distance(phy_h, method = dist_method)
    D <- as.matrix(dist_h)
    rownames(D) <- colnames(D) <- sample_names(phy_h)
    
    # meta + interaction groups
    meta_h <- as(sample_data(phy_h), "data.frame")
    meta_h$YS <- interaction(meta_h$Year, meta_h$Island_shore, sep = ".", drop = TRUE)
    meta_h <- droplevels(meta_h)
    
    groups <- levels(meta_h$YS)
    if (is.null(groups)) groups <- unique(meta_h$YS)
    if (length(groups) < 2) next
    
    combs <- utils::combn(groups, 2, simplify = FALSE)
    
    for (pair in combs) {
      g1 <- pair[1]; g2 <- pair[2]
      idx <- which(meta_h$YS %in% c(g1, g2))
      if (length(idx) < 2) next
      
      keep <- rownames(meta_h)[idx]
      M <- D[keep, keep, drop = FALSE]
      grp <- droplevels(meta_h$YS[idx])
      
      tab <- table(grp)
      if (length(tab) != 2) next
      n1 <- as.integer(tab[1]); n2 <- as.integer(tab[2])
      if (n1 < min_group_n || n2 < min_group_n) {
        out_list[[k]] <- tibble(
          Habitat = h, comparison = sprintf("%s vs %s", g1, g2),
          group1 = g1, group2 = g2,
          Year1 = strsplit(g1, "\\.")[[1]][1], Shore1 = strsplit(g1, "\\.")[[1]][2],
          Year2 = strsplit(g2, "\\.")[[1]][1], Shore2 = strsplit(g2, "\\.")[[1]][2],
          n1 = n1, n2 = n2,
          R = NA_real_, p = NA_real_, nperm_used = NA_real_
        )
        k <- k + 1
        next
      }
      
      n  <- n1 + n2
      max_perms_possible <- as.numeric(safe_choose(n, n1))
      nperm_use <- min(max_perm,
                       if (avoid_enum_msg) max(2, max_perms_possible - 1) else max_perms_possible)
      
      set.seed(seed)
      ares <- vegan::anosim(as.dist(M), grp, permutations = nperm_use)
      
      out_list[[k]] <- tibble(
        Habitat = h, comparison = sprintf("%s vs %s", g1, g2),
        group1 = g1, group2 = g2,
        Year1 = strsplit(g1, "\\.")[[1]][1], Shore1 = strsplit(g1, "\\.")[[1]][2],
        Year2 = strsplit(g2, "\\.")[[1]][1], Shore2 = strsplit(g2, "\\.")[[1]][2],
        n1 = n1, n2 = n2,
        R  = unname(ares$statistic),
        p  = unname(ares$signif),
        nperm_used = nperm_use
      )
      k <- k + 1
    }
  }
  
  if (!length(out_list)) {
    return(tibble(
      Habitat = character(), comparison = character(),
      group1 = character(), group2 = character(),
      Year1 = character(), Shore1 = character(),
      Year2 = character(), Shore2 = character(),
      n1 = integer(), n2 = integer(),
      R = double(), p = double(), nperm_used = double()
    ))
  }
  
  bind_rows(out_list) %>%
    arrange(Habitat, p, desc(R))
}
