## ------------------------------------------ ##
## Paper figure 5: LSU diversity vs drought intensification
## By: Harold Achicanoy
## WUR & ABC
## Sep 2024
## ------------------------------------------ ##

## R options and packages loading ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,MetBrewer,rnaturalearth,tidyverse,ggpubr,tseries,biscale,pals,cowplot,gridExtra,grid))

# Activate Windows fonts
extrafont::font_import()
extrafont::loadfonts(device = 'win')
extrafont::fonts()

# Functions
make_maps_legend <- function(dims = 4,
                             xlabl = 'Drought trend severity',
                             ylabl = 'Livestock units diversity',
                             xbrks = c('Negative','Low','Medium','High'),
                             ybrks = c('Very low','Low','Medium','High'),
                             palette = 'GrPink2',
                             categories = colours)
{
  categories$x <- strsplit(x = categories$bi_class, split = '-') |> purrr::map(1) |> unlist() |> as.numeric()
  categories$y <- strsplit(x = categories$bi_class, split = '-') |> purrr::map(2) |> unlist() |> as.numeric()
  legend <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = categories, mapping = ggplot2::aes(x = x, y = y, fill = bi_class)) +
    biscale::bi_scale_fill(pal = palette, dim = dims) +
    ggplot2::theme_minimal() +
    ggplot2::coord_equal() +
    ggplot2::xlab(xlabl) +
    ggplot2::ylab(ylabl) +
    ggplot2::scale_x_continuous(breaks = 1:dims, labels = xbrks) +
    ggplot2::scale_y_continuous(breaks = 1:dims, labels = ybrks) +
    ggplot2::theme(legend.position = 'none',
                   axis.text       = element_text(family = 'serif'),
                   axis.text.x     = element_text(size = 14, colour = 'black', angle = 45, hjust = 1),
                   axis.text.y     = element_text(size = 14, colour = 'black', angle = 45),
                   axis.title      = element_text(size = 16, colour = 'black', face = 'bold', family = 'serif'),
                   axis.line       = element_blank(),
                   axis.ticks      = element_blank())
  return(legend)
}

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

## Figure 5-supplementary ----
# Bivariate map of SPEI severity vs livestock units diversity (livestock exposure)

# Crop classes diversity
lsu_diversity <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/lsu_diversity_25km.tif'))

# Index characterization: extreme trend
idx_sxt <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/SPEI-6_extreme_trend.tif'))

# Zonal statistics per extreme drought cluster (mean)
agex_sgn_extreme_trend <- terra::zonal(x = idx_sxt, z = agex_sgn_num, fun = 'mean', na.rm = T, as.raster = T)
metrics1 <- terra::zonal(x = idx_sxt, z = agex_sgn_num, fun = 'mean', na.rm = T)
agex_sgn_lsu_diversity <- terra::zonal(x = lsu_diversity, z = agex_sgn_num, fun = 'mean', na.rm = T, as.raster = T)
metrics2 <- terra::zonal(x = lsu_diversity, z = agex_sgn_num, fun = 'mean', na.rm = T)
metrics <- dplyr::left_join(x = metrics1, y = metrics2, by = 'extreme_cluster'); rm(metrics1, metrics2)

# Merge data.frames
agex_sgn_severity_vs_lvst <- terra::as.data.frame(x = c(agex_sgn, agex_sgn_extreme_trend, agex_sgn_lsu_diversity), xy = T, cell = T)
rm(agex_sgn_extreme_trend, agex_sgn_lsu_diversity)

# Produce bivariate maps
# Define severity and diversity quantiles
severity_brks <- quantile(x = metrics$`SPEI-6_slope_95th`, probs = seq(0,1,1/4))
severity_brks[1] <- terra::as.data.frame(idx_sxt) |> min()
severity_brks[5] <- terra::as.data.frame(idx_sxt) |> max()
diversity_brks <- quantile(x = metrics$livestock_units_diversity, probs = seq(0,1,1/4))
diversity_brks[1] <- terra::as.data.frame(lsu_diversity)|> min()
diversity_brks[5] <- terra::as.data.frame(lsu_diversity)|> max()
# Create categories
agex_sgn_severity_vs_lvst$Severity_class <- cut(x = agex_sgn_severity_vs_lvst$`SPEI-6_slope_95th`, breaks = severity_brks) |> as.numeric() |> as.character()
agex_sgn_severity_vs_lvst$Diversity_class <- cut(x = agex_sgn_severity_vs_lvst$livestock_units_diversity, breaks = diversity_brks) |> as.numeric() |> as.character()
# Create bivariate categories
agex_sgn_severity_vs_lvst$bi_class <- paste0(agex_sgn_severity_vs_lvst$Severity_class,'-',agex_sgn_severity_vs_lvst$Diversity_class)
agex_sgn_severity_vs_lvst <- agex_sgn_severity_vs_lvst[,c('cell','x','y','extreme_cluster','SPEI-6_slope_95th','livestock_units_diversity','bi_class')]
# Define color palette for bivariate categories
colours <- data.frame(bi_class = sort(unique(agex_sgn_severity_vs_lvst$bi_class)))
colours$category <- 1:nrow(colours)
agex_sgn_severity_vs_lvst <- dplyr::left_join(x = agex_sgn_severity_vs_lvst, y = colours, by = 'bi_class')

# Create auxiliary raster of bivariate categories
aux <- agex_sgn
terra::values(aux) <- NA
aux[agex_sgn_severity_vs_lvst$cell] <- agex_sgn_severity_vs_lvst$category
names(aux) <- 'category'

# Create a list of maps (per continent)
ggm <- list()
for(i in 1:length(shp)){
  
  # Coordinates
  aux_dfm <- terra::as.data.frame(x = terra::crop(aux, terra::ext(shp[[i]])), xy = T, cell = T)
  aux_dfm <- dplyr::left_join(x = aux_dfm, y = colours, by = 'category')
  # Map
  ggm[[i]] <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = aux_dfm, aes(x, y, fill = bi_class)) +
    ggplot2::geom_sf(data = sf::st_as_sf(shp[[i]]), fill = NA) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = 'none', legend.key.width = unit(1.3, 'cm')) +
    biscale::bi_scale_fill(pal = 'GrPink2', dim = 4)
  
}; rm(i,aux_dfm)

### Legend ----
fig3_leg <- make_maps_legend(dims = 4,
                             xlabl = 'Extreme drought\nintensification',
                             ylabl = 'Livestock units\ndiversity',
                             xbrks = c('Negative','Low','Medium','High'),
                             ybrks = c('Very low','Low','Medium','High'),
                             palette = 'GrPink2',
                             categories = colours) # ggplot2::ggsave(filename = paste0('D:/Figure3_paper1_legend.png'), plot = fig3_leg, device = 'png', width = 5, height = 5, units = 'in', dpi = 350)

### Panel a) North America ----
# Cluster: 227 (USA, Mexico)
# One growing season

idx_ts <- terra::rast(paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_25km/one_s1_spei-6_25km.tif'))
# Filter time series per cluster of interest
aux <- agex_sgn; aux[aux != 227] <- NA
idx_ts_msk <- terra::mask(x = idx_ts, mask = aux) |> terra::trim()
# Time series plotting
idx_ts_msk_dfm <- terra::as.data.frame(x = idx_ts_msk, xy = T, cell = F)
gg_ts_NA <- idx_ts_msk_dfm |>
  tidyr::pivot_longer(cols = 3:ncol(idx_ts_msk_dfm), names_to = 'Years', values_to = 'SPEI-6') |>
  dplyr::mutate(Years = as.numeric(gsub('spei-6_', '', Years)),
                ID = paste0(x, '-', y)) |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = Years, y = `SPEI-6`, group = ID)) +
  ggplot2::geom_line(alpha = 0.1) +
  ggplot2::geom_hline(yintercept = 1.5, linetype = 'dashed', color = 'red', size = 1.1) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'plain', family = 'serif'))

# Add polygon extent to the corresponding map
ext_vct_NA <- as.vector(terra::ext(idx_ts_msk))
ggm_NA <- ggm[[1]] +
  ggplot2::geom_rect(aes(xmin = ext_vct_NA[1],
                         xmax = ext_vct_NA[2],
                         ymin = ext_vct_NA[3],
                         ymax = ext_vct_NA[4]),
                     color = 'red', size = 1.5, fill = NA)

### Panel b) South America ----
# Cluster: 402 (Bolivia,Peru,Argentina)
# One growing season

idx_ts <- terra::rast(paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_25km/one_s1_spei-6_25km.tif'))
# Filter time series per cluster of interest
aux <- agex_sgn; aux[aux != 402] <- NA
idx_ts_msk <- terra::mask(x = idx_ts, mask = aux) |> terra::trim()
# Time series plotting
idx_ts_msk_dfm <- terra::as.data.frame(x = idx_ts_msk, xy = T, cell = F)
gg_ts_SA <- idx_ts_msk_dfm |>
  tidyr::pivot_longer(cols = 3:ncol(idx_ts_msk_dfm), names_to = 'Years', values_to = 'SPEI-6') |>
  dplyr::mutate(Years = as.numeric(gsub('spei-6_', '', Years)),
                ID = paste0(x, '-', y)) |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = Years, y = `SPEI-6`, group = ID)) +
  ggplot2::geom_line(alpha = 0.1) +
  ggplot2::geom_hline(yintercept = 1.5, linetype = 'dashed', color = 'red', size = 1.1) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'plain', family = 'serif'))

# Add polygon extent to the corresponding map
ext_vct_SA <- as.vector(terra::ext(idx_ts_msk))
ggm_SA <- ggm[[2]] +
  ggplot2::geom_rect(aes(xmin = ext_vct_SA[1],
                         xmax = ext_vct_SA[2],
                         ymin = ext_vct_SA[3],
                         ymax = ext_vct_SA[4]),
                     color = 'red', size = 1.5, fill = NA)

### Panel c) Africa ----
# Cluster: 289 (Nigeria,Benin,Togo,Ghana)
# Two growing season. Chosen season S1

idx_ts <- terra::rast(paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_25km/two_s1_spei-6_25km.tif'))
# Filter time series per cluster of interest
aux <- agex_sgn; aux[aux != 289] <- NA
idx_ts_msk <- terra::mask(x = idx_ts, mask = aux) |> terra::trim()
# Time series plotting
idx_ts_msk_dfm <- terra::as.data.frame(x = idx_ts_msk, xy = T, cell = F)
gg_ts_AF <- idx_ts_msk_dfm |>
  tidyr::pivot_longer(cols = 3:ncol(idx_ts_msk_dfm), names_to = 'Years', values_to = 'SPEI-6') |>
  dplyr::mutate(Years = as.numeric(gsub('spei-6_', '', Years)),
                ID = paste0(x, '-', y)) |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = Years, y = `SPEI-6`, group = ID)) +
  ggplot2::geom_line(alpha = 0.1) +
  ggplot2::geom_hline(yintercept = 1.5, linetype = 'dashed', color = 'red', size = 1.1) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'plain', family = 'serif'))

# Add polygon extent to the corresponding map
ext_vct_AF <- as.vector(terra::ext(idx_ts_msk))
ggm_AF <- ggm[[3]] +
  ggplot2::geom_rect(aes(xmin = ext_vct_AF[1],
                         xmax = ext_vct_AF[2],
                         ymin = ext_vct_AF[3],
                         ymax = ext_vct_AF[4]),
                     color = 'red', size = 1.5, fill = NA)

### Panel d) Europe ----
# Cluster: 50 (Ukraine,Russia)
# One growing season

idx_ts <- terra::rast(paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_25km/one_s1_spei-6_25km.tif'))
# Filter time series per cluster of interest
aux <- agex_sgn; aux[aux != 50] <- NA
idx_ts_msk <- terra::mask(x = idx_ts, mask = aux) |> terra::trim()
# Time series plotting
idx_ts_msk_dfm <- terra::as.data.frame(x = idx_ts_msk, xy = T, cell = F)
gg_ts_EU <- idx_ts_msk_dfm |>
  tidyr::pivot_longer(cols = 3:ncol(idx_ts_msk_dfm), names_to = 'Years', values_to = 'SPEI-6') |>
  dplyr::mutate(Years = as.numeric(gsub('spei-6_', '', Years)),
                ID = paste0(x, '-', y)) |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = Years, y = `SPEI-6`, group = ID)) +
  ggplot2::geom_line(alpha = 0.1) +
  ggplot2::geom_hline(yintercept = 1.5, linetype = 'dashed', color = 'red', size = 1.1) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'plain', family = 'serif'))

# Add polygon extent to the corresponding map
ext_vct_EU <- as.vector(terra::ext(idx_ts_msk))
ggm_EU <- ggm[[4]] +
  ggplot2::geom_rect(aes(xmin = ext_vct_EU[1],
                         xmax = ext_vct_EU[2],
                         ymin = ext_vct_EU[3],
                         ymax = ext_vct_EU[4]),
                     color = 'red', size = 1.5, fill = NA)

### Panel e) Asia ----
# Cluster: 171 (People's Republic of China,South Korea,Japan)
# One growing season

idx_ts <- terra::rast(paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_25km/one_s1_spei-6_25km.tif'))
# Filter time series per cluster of interest
aux <- agex_sgn; aux[aux != 171] <- NA
idx_ts_msk <- terra::mask(x = idx_ts, mask = aux) |> terra::trim()
# Time series plotting
idx_ts_msk_dfm <- terra::as.data.frame(x = idx_ts_msk, xy = T, cell = F)
gg_ts_AS <- idx_ts_msk_dfm |>
  tidyr::pivot_longer(cols = 3:ncol(idx_ts_msk_dfm), names_to = 'Years', values_to = 'SPEI-6') |>
  dplyr::mutate(Years = as.numeric(gsub('spei-6_', '', Years)),
                ID = paste0(x, '-', y)) |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = Years, y = `SPEI-6`, group = ID)) +
  ggplot2::geom_line(alpha = 0.1) +
  ggplot2::geom_hline(yintercept = 1.5, linetype = 'dashed', color = 'red', size = 1.1) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'plain', family = 'serif'))

# Add polygon extent to the corresponding map
ext_vct_AS <- as.vector(terra::ext(idx_ts_msk))
ggm_AS <- ggm[[5]] +
  ggplot2::geom_rect(aes(xmin = ext_vct_AS[1],
                         xmax = ext_vct_AS[2],
                         ymin = ext_vct_AS[3],
                         ymax = ext_vct_AS[4]),
                     color = 'red', size = 1.5, fill = NA)

### Panel f) Oceania ----
# Cluster: 444 (Australia)
# One growing season

idx_ts <- terra::rast(paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_25km/one_s1_spei-6_25km.tif'))
# Filter time series per cluster of interest
aux <- agex_sgn; aux[aux != 444] <- NA
idx_ts_msk <- terra::mask(x = idx_ts, mask = aux) |> terra::trim()
# Time series plotting
idx_ts_msk_dfm <- terra::as.data.frame(x = idx_ts_msk, xy = T, cell = F)
gg_ts_OC <- idx_ts_msk_dfm |>
  tidyr::pivot_longer(cols = 3:ncol(idx_ts_msk_dfm), names_to = 'Years', values_to = 'SPEI-6') |>
  dplyr::mutate(Years = as.numeric(gsub('spei-6_', '', Years)),
                ID = paste0(x, '-', y)) |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = Years, y = `SPEI-6`, group = ID)) +
  ggplot2::geom_line(alpha = 0.1) +
  ggplot2::geom_hline(yintercept = 1.5, linetype = 'dashed', color = 'red', size = 1.1) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'plain', family = 'serif'))

# Add polygon extent to the corresponding map
ext_vct_OC <- as.vector(terra::ext(idx_ts_msk))
ggm_OC <- ggm[[6]] +
  ggplot2::geom_rect(aes(xmin = ext_vct_OC[1],
                         xmax = ext_vct_OC[2],
                         ymin = ext_vct_OC[3],
                         ymax = ext_vct_OC[4]),
                     color = 'red', size = 1.5, fill = NA)

## Full figure ----
layout <- matrix(c(1:12,NA,13,13,NA), nrow = 4, ncol = 4, byrow = T)

sts <- agex_sgn_metrics[agex_sgn_metrics$extreme_cluster %in% c(157,420,270,33,201,444),]
sts[,c('harvested_area_total','harvested_area_min','harvested_area_max')]

fg_01 <- ggpubr::annotate_figure(gg_ts_NA,
                                 top = text_grob(
                                   label  = 'a) Transnational. One-GS',
                                   face   = 'plain',
                                   size   = 13,
                                   family = 'serif'
                                 ))
fg_02 <- ggpubr::annotate_figure(ggm_NA,
                                 top = text_grob(
                                   label  = 'a) North America',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_03 <- ggpubr::annotate_figure(ggm_SA,
                                 top = text_grob(
                                   label  = 'b) South America',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_04 <- ggpubr::annotate_figure(gg_ts_SA,
                                 top = text_grob(
                                   label  = 'b) Transnational. One-GS',
                                   face   = 'plain',
                                   size   = 13,
                                   family = 'serif'
                                 ))
fg_05 <- ggpubr::annotate_figure(gg_ts_AF,
                                 top = text_grob(
                                   label  = 'c) Transnational. Two-GS',
                                   face   = 'plain',
                                   size   = 13,
                                   family = 'serif'
                                 ))
fg_06 <- ggpubr::annotate_figure(ggm_AF,
                                 top = text_grob(
                                   label  = 'c) Africa',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_07 <- ggpubr::annotate_figure(ggm_EU,
                                 top = text_grob(
                                   label  = 'd) Europe',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_08 <- ggpubr::annotate_figure(gg_ts_EU,
                                 top = text_grob(
                                   label  = 'd) Transnational. One-GS',
                                   face   = 'plain',
                                   size   = 13,
                                   family = 'serif'
                                 ))
fg_09 <- ggpubr::annotate_figure(gg_ts_AS,
                                 top = text_grob(
                                   label  = 'e) Transnational. One-GS',
                                   face   = 'plain',
                                   size   = 13,
                                   family = 'serif'
                                 ))
fg_10 <- ggpubr::annotate_figure(ggm_AS,
                                 top = text_grob(
                                   label  = 'e) Asia',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_11 <- ggpubr::annotate_figure(ggm_OC,
                                 top = text_grob(
                                   label  = 'f) Oceania',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_12 <- ggpubr::annotate_figure(gg_ts_OC,
                                 top = text_grob(
                                   label  = 'f) National. One-GS',
                                   face   = 'plain',
                                   size   = 13,
                                   family = 'serif'
                                 ))
fg_13 <- ggpubr::annotate_figure(fig3_leg, top = NULL)

fig3 <- gridExtra::grid.arrange(fg_01, fg_02, fg_03, fg_04,
                                fg_05, fg_06, fg_07, fg_08,
                                fg_09, fg_10, fg_11, fg_12,
                                fg_13,
                                layout_matrix = layout)
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/Figure5_paper1-livestock_diversity.png'), plot = fig3, device = 'png', width = 12, height = 10, units = 'in', dpi = 350)
rm(ggm,
   ggm_NA,ggm_SA,ggm_AF,ggm_EU,ggm_AS,ggm_OC,
   gg_ts_NA,gg_ts_SA,gg_ts_AF,gg_ts_EU,gg_ts_AS,gg_ts_OC,
   fg_01,fg_02,fg_03,fg_04,fg_05,fg_06,fg_07,fg_08,fg_09,fg_10,fg_11,fg_12,fg_13,
   fig3_leg,fig3)
