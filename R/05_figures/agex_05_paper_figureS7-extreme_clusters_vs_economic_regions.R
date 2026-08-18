# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Figure S7: Extreme agricultural drought clusters overlaid with
# economic comunity regions
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
agex_sgn <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn) <- 'extreme_cluster'
cls <- data.frame(extreme_cluster = sort(unique(terra::values(agex_sgn))))
cls$id <- cls$extreme_cluster
cls <- cls[,c('id','extreme_cluster')]
cls$extreme_cluster <- as.character(cls$extreme_cluster)
levels(agex_sgn) <- cls; rm(cls)

## Shapefiles ----
# Shapefiles ----
sam <- rnaturalearth::ne_countries(scale = 'large', continent = 'south america', returnclass = 'sv')
sam <- terra::aggregate(sam); sam$Continent <- 'South America'
esam <- terra::vect("D:/hotspots_agriculture/outputs/rec_gpkg/MERCOSUR_Southern_Common_Market.gpkg")
esam <- terra::aggregate(esam); esam$Community <- 'MERCOSUR'
sam <- terra::crop(x = sam, y = terra::ext(esam))
asi <- rnaturalearth::ne_countries(scale = 'large', continent = 'asia', returnclass = 'sv')
asi <- terra::aggregate(asi); asi$Continent <- 'Asia'
easi <- terra::vect("D:/hotspots_agriculture/outputs/rec_gpkg/ASEAN_Association_of_Southeast_Asian_Nations.gpkg")
easi <- terra::aggregate(easi); easi$Community <- 'ASEAN'
asi <- terra::crop(x = asi, y = terra::ext(easi))
eur <- rnaturalearth::ne_countries(scale = 'large', continent = 'europe', returnclass = 'sv')
eur <- terra::aggregate(eur); eur$Continent <- 'Europe'
eeur <- terra::vect("D:/hotspots_agriculture/outputs/rec_gpkg/EU_European_Union.gpkg")
eeur <- terra::aggregate(eeur); eeur$Community <- 'EU'
eur <- terra::crop(x = eur, y = terra::ext(eeur))
afr <- rnaturalearth::ne_countries(scale = 'large', continent = 'africa', returnclass = 'sv')
afr <- terra::aggregate(afr); afr$Continent <- 'Africa'
eafr <- terra::vect("D:/hotspots_agriculture/outputs/rec_gpkg/SADC_Southern_African_Development_Community.gpkg")
eafr <- terra::aggregate(eafr); eafr$Community <- 'SADC'
afr <- terra::crop(x = afr, y = terra::ext(eafr))
shp <- rbind(sam, asi, eur, afr)
rm(sam, asi, eur, afr)
region <- rbind(esam, easi, eeur, eafr)
rm(esam, easi, eeur, eafr)

csam <- terra::vect("D:/hotspots_agriculture/outputs/rec_gpkg/MERCOSUR_Southern_Common_Market.gpkg")
casi <- terra::vect("D:/hotspots_agriculture/outputs/rec_gpkg/ASEAN_Association_of_Southeast_Asian_Nations.gpkg")
ceur <- terra::vect("D:/hotspots_agriculture/outputs/rec_gpkg/EU_European_Union.gpkg")
cafr <- terra::vect("D:/hotspots_agriculture/outputs/rec_gpkg/SADC_Southern_African_Development_Community.gpkg")
countries <- list(csam, casi, ceur, cafr)

# Figure S7 ----
n_col    <- length(unique(terra::values(agex_sgn,na.rm = T))) # Number of unique colors needed
col_pltt <- MetBrewer::met.brewer(name = 'Signac', n = n_col) # Color palette
set.seed(1235); col_pltt <- sample(x = col_pltt, size = n_col, replace = F); rm(n_col) # Randomize colors assignment

# Create a list of maps (per continent)
ggm <- list()
for(i in 1:length(shp)){
  
  # Coordinates
  agex_sgn_dfm <- terra::as.data.frame(x = terra::crop(agex_sgn, terra::ext(shp[i])), xy = T, cell = T)
  # Setting-up colors
  stp <- data.frame(extreme_cluster = factor(sort(unique(terra::values(agex_sgn, na.rm = T)))), color = col_pltt)
  agex_sgn_dfm <- dplyr::left_join(x = agex_sgn_dfm, y = stp, by = 'extreme_cluster')
  agex_sgn_dfm$extreme_cluster <- as.factor(agex_sgn_dfm$extreme_cluster)
  sub_stp <- unique(agex_sgn_dfm[,c('extreme_cluster','color')])
  # Map
  ggm[[i]] <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = agex_sgn_dfm, aes(x, y, fill = extreme_cluster)) +
    ggplot2::geom_sf(data = sf::st_as_sf(region[i]), fill = 'orange', alpha = 0.5, color = 'black', linewidth = 0.3) +
    ggplot2::geom_sf(data = sf::st_as_sf(shp[i]), fill = NA) +
    ggplot2::geom_sf(data = sf::st_as_sf(countries[[i]]), fill = NA, color = 'black', linewidth = 0.1) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = 'none', legend.key.width = unit(1.3, 'cm')) +
    ggplot2::scale_fill_manual(values = sub_stp$color, breaks = sub_stp$extreme_cluster)
  
}; rm(i,stp,sub_stp,n_col,col_pltt)

# Organize maps into one graph adding panel labels (alternative 2)
fg_01 <- ggpubr::annotate_figure(ggm[[1]], top = text_grob(label = 'a) MERCOSUR', face = 'plain', size = 20, family = 'serif'))
fg_02 <- ggpubr::annotate_figure(ggm[[2]], top = text_grob(label = 'b) ASEAN', face = 'plain', size = 20, family = 'serif'))
fg_03 <- ggpubr::annotate_figure(ggm[[3]], top = text_grob(label = 'c) EU', face = 'plain', size = 20, family = 'serif'))
fg_04 <- ggpubr::annotate_figure(ggm[[4]], top = text_grob(label = 'd) SADC', face = 'plain', size = 20, family = 'serif'))

layout <- matrix(1:4, nrow = 2, ncol = 2, byrow = T)

figS7 <- gridExtra::grid.arrange(fg_01, fg_02,
                                 fg_03, fg_04,
                                 layout_matrix = layout); rm(ggm,fg_01,fg_02,fg_03,fg_04,layout)
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/FigureS7-agex-extreme_drought_clusters-economic_regions.png'), plot = figS7, device = 'png', width = 10, height = 10, units = 'in', dpi = 350); rm(figS7)
