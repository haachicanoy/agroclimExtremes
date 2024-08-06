# ------------------------------------------ #
# Agro-ecological zones plotting
# By: Harold Achicanoy
# WUR & ABC
# Aug 2024
# ------------------------------------------ #

## R options and packages loading ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,rnaturalearth,tidyverse))

root <- '//CATALOGUE/WFP_ClimateRiskPr1/agroclimExtremes'

## World shapefile ----
wrld <- rnaturalearth::ne_coastline(scale = 'medium', returnclass = 'sv')

## Extreme drought clusters ----
# Shapefile
agex <- terra::vect(paste0(root,'/agex_results/agex_results_clusters/vct_agex_global_spei-6.gpkg'))
# Categorical
agex_sgn <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn) <- 'extreme_cluster'
cls <- data.frame(extreme_cluster = sort(unique(terra::values(agex_sgn))))
cls$id <- cls$extreme_cluster
cls <- cls[,c('id','extreme_cluster')]
cls$extreme_cluster <- as.character(cls$extreme_cluster)
levels(agex_sgn) <- cls; rm(cls)

# Numerical
agex_sgn_num <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn_num) <- 'extreme_cluster'

# Color palette definition
n_col    <- length(unique(terra::values(agex_sgn,na.rm = T))) # Number of unique colors needed
col_pltt <- MetBrewer::met.brewer(name = 'Signac', n = n_col) # Color palette
set.seed(1235); col_pltt <- sample(x = col_pltt, size = n_col, replace = F); rm(n_col) # Randomize colors assignment

# Coordinates
agex_sgn_dfm <- terra::as.data.frame(x = agex_sgn, xy = T, cell = T)
# Setting-up colors
stp <- data.frame(extreme_cluster = factor(sort(unique(terra::values(agex_sgn, na.rm = T)))), color = col_pltt)
agex_sgn_dfm <- dplyr::left_join(x = agex_sgn_dfm, y = stp, by = 'extreme_cluster')
agex_sgn_dfm$extreme_cluster <- as.factor(agex_sgn_dfm$extreme_cluster)
sub_stp <- unique(agex_sgn_dfm[,c('extreme_cluster','color')])
# Map
ggm <- ggplot2::ggplot() +
  ggplot2::geom_tile(data = agex_sgn_dfm, aes(x, y, fill = extreme_cluster)) +
  ggplot2::geom_sf(data = sf::st_as_sf(wrld), fill = NA) +
  ggplot2::coord_sf() +
  # ggplot2::theme_void() +
  # ggplot2::theme(legend.position = 'none', legend.key.width = unit(1.3, 'cm')) +
  ggplot2::scale_fill_manual(values = sub_stp$color, breaks = sub_stp$extreme_cluster) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = 'none')
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/supplementary/agex.png'), plot = ggm, width = 10, height = 5, units = 'in', dpi = 350)

## Agro-ecological zones ----
gaez <- terra::rast(paste0(root,'/agex_raw_data/gaez_v4/aez/aez_v9v2_CRUTS32_Hist_8110_100_avg.tif'))
gaez_masked <- terra::mask(x = gaez, mask = agex)

## Save results ----
png(filename = paste0(root,'/agex_results/agex_figures/supplementary/gaez.png'), width = 10, height = 5, units = 'in', res = 350)
plot(wrld)
plot(gaez, add = T)
dev.off()

png(filename = paste0(root,'/agex_results/agex_figures/supplementary/gaez_mskd.png'), width = 10, height = 5, units = 'in', res = 350)
plot(wrld)
plot(gaez_masked, add = T)
dev.off()

png(filename = paste0(root,'/agex_results/agex_figures/supplementary/agex.png'), width = 10, height = 5, units = 'in', res = 350)
plot(wrld)
plot(agex_sgn_num, col = stp$color, add = T, legend = F)
dev.off()
