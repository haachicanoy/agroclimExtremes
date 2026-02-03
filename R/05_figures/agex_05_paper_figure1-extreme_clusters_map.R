# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Figure 1: Extreme agricultural drought clusters
# By: Harold Achicanoy
# WUR & ABC
# Created in October 2024
# Modified in February 2026
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

# Metrics
agex_sgn_metrics <- utils::read.csv(paste0(root,'/agex_results/agex_all_metrics.csv'))

## Shapefiles ----
wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
wrl <- wrl[,c('adm0_iso','name_en','continent','subregion')]

# Continent shapefiles. Some of them cropped
afr <- rnaturalearth::ne_countries(scale = 'large', continent = 'africa', returnclass = 'sv')
eur <- rnaturalearth::ne_countries(scale = 'large', continent = 'europe', returnclass = 'sv')
eur <- terra::crop(x = eur, y = terra::ext(c(-20,100,30,75))) # xmin, xmax = 150, ymin, ymax
asi <- rnaturalearth::ne_countries(scale = 'large', continent = 'asia', returnclass = 'sv')
oce <- rnaturalearth::ne_countries(scale = 'large', continent = 'oceania', returnclass = 'sv')
oce <- terra::crop(x = oce, y = terra::ext(c(90,180,-54.75,20.55)))
nam <- rnaturalearth::ne_countries(scale = 'large', continent = 'north america', returnclass = 'sv')
nam <- terra::crop(x = nam, y = terra::ext(c(-140,-50,-0.39,60)))
sam <- rnaturalearth::ne_countries(scale = 'large', continent = 'south america', returnclass = 'sv')

# Put them all together
shp <- list(nam, sam, afr, eur, asi, oce); rm(nam, sam, afr, eur, asi, oce)

# Figure 1 ----
# Extreme drought clusters per continent
# Color palette definition
n_col    <- length(unique(terra::values(agex_sgn,na.rm = T))) # Number of unique colors needed
col_pltt <- MetBrewer::met.brewer(name = 'Signac', n = n_col) # Color palette
set.seed(1235); col_pltt <- sample(x = col_pltt, size = n_col, replace = F); rm(n_col) # Randomize colors assignment

# Create a list of maps (per continent)
ggm <- list()
for(i in 1:length(shp)){
  
  # Coordinates
  agex_sgn_dfm <- terra::as.data.frame(x = terra::crop(agex_sgn, terra::ext(shp[[i]])), xy = T, cell = T)
  # Setting-up colors
  stp <- data.frame(extreme_cluster = factor(sort(unique(terra::values(agex_sgn, na.rm = T)))), color = col_pltt)
  agex_sgn_dfm <- dplyr::left_join(x = agex_sgn_dfm, y = stp, by = 'extreme_cluster')
  agex_sgn_dfm$extreme_cluster <- as.factor(agex_sgn_dfm$extreme_cluster)
  sub_stp <- unique(agex_sgn_dfm[,c('extreme_cluster','color')])
  # Map
  ggm[[i]] <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = agex_sgn_dfm, aes(x, y, fill = extreme_cluster)) +
    ggplot2::geom_sf(data = sf::st_as_sf(shp[[i]]), fill = NA) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = 'none', legend.key.width = unit(1.3, 'cm')) +
    ggplot2::scale_fill_manual(values = sub_stp$color, breaks = sub_stp$extreme_cluster)
  
}; rm(i,stp,sub_stp,n_col,col_pltt)

# Organize maps into one graph adding panel labels (alternative 2)
fg_01 <- ggpubr::annotate_figure(ggm[[1]], top = text_grob(label = 'a) North America', face = 'plain', size = 20, family = 'serif'))
fg_02 <- ggpubr::annotate_figure(ggm[[2]], top = text_grob(label = 'b) South America', face = 'plain', size = 20, family = 'serif'))
fg_03 <- ggpubr::annotate_figure(ggm[[3]], top = text_grob(label = 'c) Africa', face = 'plain', size = 20, family = 'serif'))
fg_04 <- ggpubr::annotate_figure(ggm[[4]], top = text_grob(label = 'd) Europe', face = 'plain', size = 20, family = 'serif'))
fg_05 <- ggpubr::annotate_figure(ggm[[5]], top = text_grob(label = 'e) Asia', face = 'plain', size = 20, family = 'serif'))
fg_06 <- ggpubr::annotate_figure(ggm[[6]], top = text_grob(label = 'f) Oceania', face = 'plain', size = 20, family = 'serif'))

layout <- matrix(1:6, nrow = 3, ncol = 2, byrow = T)

fig1 <- gridExtra::grid.arrange(fg_01, fg_02,
                                fg_03, fg_04,
                                fg_05, fg_06,
                                layout_matrix = layout); rm(ggm,fg_01,fg_02,fg_03,fg_04,fg_05,fg_06,layout)
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/Figure1-agex-extreme_drought_clusters.png'), plot = fig1, device = 'png', width = 10, height = 12.5, units = 'in', dpi = 350); rm(fig1)
