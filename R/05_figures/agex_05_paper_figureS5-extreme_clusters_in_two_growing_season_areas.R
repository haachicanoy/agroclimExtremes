# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Figure S5: Extreme agricultural drought clusters in areas with 2
# growing seasons
# By: Harold Achicanoy
# WUR & ABC
# Created in August 2026
# --------------------------------------------------------------- #

# R options and user-defined functions ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,MetBrewer,rnaturalearth,tidyverse,ggpubr,tseries,biscale,pals,cowplot,gridExtra,grid))

# Activate Windows fonts
extrafont::font_import()
extrafont::loadfonts(device = 'win')
extrafont::fonts()

## Setup arguments ----
root   <- '//CATALOGUE/AgroclimExtremes'
index  <- 'spei-6'

## Extreme drought clusters ----
# Categorical
agex_sgn1 <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_two_s1_fmadogram_k225.tif'))
names(agex_sgn1) <- 'extreme_cluster'
cls <- data.frame(extreme_cluster = sort(unique(terra::values(agex_sgn1))))
cls$id <- cls$extreme_cluster
cls <- cls[,c('id','extreme_cluster')]
cls$extreme_cluster <- as.character(cls$extreme_cluster)
levels(agex_sgn1) <- cls; rm(cls)
agex_sgn2 <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_two_s2_fmadogram_k222.tif'))
names(agex_sgn2) <- 'extreme_cluster'
cls <- data.frame(extreme_cluster = sort(unique(terra::values(agex_sgn2))))
cls$id <- cls$extreme_cluster
cls <- cls[,c('id','extreme_cluster')]
cls$extreme_cluster <- as.character(cls$extreme_cluster)
levels(agex_sgn2) <- cls; rm(cls)

# # Numerical
# agex_sgn_num <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
# names(agex_sgn_num) <- 'extreme_cluster'
# 
# # Metrics
# agex_sgn_metrics <- utils::read.csv(paste0(root,'/agex_results/agex_all_metrics.csv'))

## Shapefiles ----
# Shapefiles ----
afr <- rnaturalearth::ne_countries(scale = 'large', continent = 'africa', returnclass = 'sv')
eur <- rnaturalearth::ne_countries(scale = 'large', continent = 'europe', returnclass = 'sv')
asi <- rnaturalearth::ne_countries(scale = 'large', continent = 'asia', returnclass = 'sv')
oce <- rnaturalearth::ne_countries(scale = 'large', continent = 'oceania', returnclass = 'sv')
nam <- rnaturalearth::ne_countries(scale = 'large', continent = 'north america', returnclass = 'sv')
sam <- rnaturalearth::ne_countries(scale = 'large', continent = 'south america', returnclass = 'sv')
afr <- terra::aggregate(afr); afr$Continent <- 'Africa'
eur <- terra::aggregate(eur); eur$Continent <- 'Europe'
asi <- terra::aggregate(asi); asi$Continent <- 'Asia'
oce <- terra::aggregate(oce); oce$Continent <- 'Oceania'
nam <- terra::aggregate(nam); nam$Continent <- 'North America'
sam <- terra::aggregate(sam); sam$Continent <- 'South America'
shp <- rbind(afr, eur, asi, oce, nam, sam)
rm(afr, eur, asi, oce, nam, sam)

# Figure S5 ----
# Extreme drought clusters per growing season
# Color palette definition
n_col1    <- length(unique(terra::values(agex_sgn1, na.rm = T))) # Number of unique colors needed
n_col2    <- length(unique(terra::values(agex_sgn2, na.rm = T))) # Number of unique colors needed
col_pltt1 <- MetBrewer::met.brewer(name = 'Signac', n = n_col1) # Color palette
col_pltt2 <- MetBrewer::met.brewer(name = 'Signac', n = n_col2) # Color palette
set.seed(1235); col_pltt1 <- sample(x = col_pltt1, size = n_col1, replace = F); rm(n_col1) # Randomize colors assignment
set.seed(1234); col_pltt2 <- sample(x = col_pltt2, size = n_col2, replace = F); rm(n_col2) # Randomize colors assignment

# Create a list of maps (per continent)
ggm <- list()
# Coordinates
agex_sgn1_dfm <- terra::as.data.frame(x = agex_sgn1, xy = T, cell = T)
agex_sgn2_dfm <- terra::as.data.frame(x = agex_sgn2, xy = T, cell = T)

# Setting-up colors
stp1 <- data.frame(extreme_cluster = factor(sort(unique(terra::values(agex_sgn1, na.rm = T)))), color = col_pltt1)
agex_sgn1_dfm <- dplyr::left_join(x = agex_sgn1_dfm, y = stp1, by = 'extreme_cluster')
agex_sgn1_dfm$extreme_cluster <- as.factor(agex_sgn1_dfm$extreme_cluster)
sub_stp1 <- unique(agex_sgn1_dfm[,c('extreme_cluster','color')])
# Map
ggm[[1]] <- ggplot2::ggplot() +
  ggplot2::geom_tile(data = agex_sgn1_dfm, aes(x, y, fill = extreme_cluster)) +
  ggplot2::geom_sf(data = sf::st_as_sf(shp), fill = NA) +
  ggplot2::coord_sf() +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = 'none', legend.key.width = unit(1.3, 'cm')) +
  ggplot2::scale_fill_manual(values = sub_stp1$color, breaks = sub_stp1$extreme_cluster)

# Setting-up colors
stp2 <- data.frame(extreme_cluster = factor(sort(unique(terra::values(agex_sgn2, na.rm = T)))), color = col_pltt2)
agex_sgn2_dfm <- dplyr::left_join(x = agex_sgn2_dfm, y = stp2, by = 'extreme_cluster')
agex_sgn2_dfm$extreme_cluster <- as.factor(agex_sgn2_dfm$extreme_cluster)
sub_stp2 <- unique(agex_sgn2_dfm[,c('extreme_cluster','color')])
# Map
ggm[[2]] <- ggplot2::ggplot() +
  ggplot2::geom_tile(data = agex_sgn2_dfm, aes(x, y, fill = extreme_cluster)) +
  ggplot2::geom_sf(data = sf::st_as_sf(shp), fill = NA) +
  ggplot2::coord_sf() +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = 'none', legend.key.width = unit(1.3, 'cm')) +
  ggplot2::scale_fill_manual(values = sub_stp2$color, breaks = sub_stp2$extreme_cluster)

# Organize maps into one graph adding panel labels (alternative 2)
fg_01 <- ggpubr::annotate_figure(ggm[[1]], top = text_grob(label = 'Growing season 1', face = 'plain', size = 15, family = 'serif'))
fg_02 <- ggpubr::annotate_figure(ggm[[2]], top = text_grob(label = 'Growing season 2', face = 'plain', size = 15, family = 'serif'))

layout <- matrix(1:2, nrow = 2, ncol = 1, byrow = T)

fig1 <- gridExtra::grid.arrange(fg_01,
                                fg_02, layout_matrix = layout); rm(ggm,fg_01,fg_02,layout)
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/FigureS5-agex-extreme_drought_clusters_gs2_areas.png'), plot = fig1, device = 'png', width = 10, height = 9, units = 'in', dpi = 350); rm(fig1)
