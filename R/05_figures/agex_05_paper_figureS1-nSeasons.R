# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Figure S1: number of growing seasons map
# By: Harold Achicanoy
# WUR & ABC
# Created in 2025
# Modified in February 2026
# --------------------------------------------------------------- #

# R options and user-defined functions ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,sf,rnaturalearth,tidyverse,ggpubr,tseries,biscale,pals,cowplot,gridExtra,grid))

# Setup arguments ----
root   <- '//CATALOGUE/AgroclimExtremes'
index  <- 'spei-6'

# Growing season places ----
nseasons <- terra::rast(file.path(root,'agex_raw_data/agex_nseasons_25km_corrected.tif'))

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
wrl <- rbind(afr, eur, asi, oce, nam, sam)
rm(afr, eur, asi, oce, nam, sam)

# Geographical coordinates ----
aux_dfm <- terra::as.data.frame(x = nseasons, xy = T, cell = T)
names(aux_dfm)[ncol(aux_dfm)] <- 'nSeasons'
aux_dfm$nSeasons <- factor(aux_dfm$nSeasons)

# Map ----
ggm <- ggplot2::ggplot() +
  ggplot2::geom_tile(data = aux_dfm, aes(x, y, fill = nSeasons)) +
  ggplot2::geom_sf(data = sf::st_as_sf(wrl), fill = NA) +
  ggplot2::coord_sf() +
  ggplot2::theme_void() +
  ggplot2::labs(fill = 'Number of growing seasons') +
  ggplot2::theme(legend.title    = element_text(family = 'serif', face = 'bold', colour = 'black', size = 17),
                 legend.text     = element_text(family = 'serif', colour = 'black', size = 15),
                 legend.position = 'bottom', legend.key.width = unit(1.3, 'cm')); ggm
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/FigureS1-agex-nSeasons.png'), plot = ggm, device = 'png', width = 12.5, height = 10, units = 'in', dpi = 350)
