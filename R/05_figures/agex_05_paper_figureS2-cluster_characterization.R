# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Figure S2: clusters characterization
# By: Harold Achicanoy
# WUR & ABC
# Created in 2025
# Modified in February 2026
# --------------------------------------------------------------- #

options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(terra,MetBrewer,sf,rnaturalearth,tidyverse,ggpubr,tseries,biscale,pals,cowplot,gridExtra,grid))

# Setup arguments ----
root   <- '//CATALOGUE/AgroclimExtremes'
index  <- 'spei-6'

# Extreme drought clusters ----
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

# Shapefiles ----
# Continent shapefiles
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

# Number of countries spanning clusters ----
# Number of countries per cluster
countries_per_cluster <- agex_sgn_num
countries_per_cluster <- terra::classify(x = countries_per_cluster, as.matrix(agex_sgn_metrics[,c('extreme_cluster','countries_count')]))
ncountries <- terra::zonal(x = countries_per_cluster, z = wrl, fun = 'mean', na.rm = T); rm(countries_per_cluster)
ncountries$Continent <- wrl$Continent
names(ncountries)[1] <- 'Average'
gg1 <- ncountries |>
  ggplot2::ggplot(aes(x = reorder(Continent, Average), y = Average)) +
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::xlab('') +
  ggplot2::ylab('Number of countries by cluster') +
  ggplot2::coord_flip() +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'plain', family = 'serif'))

# Drought intensification ----
# Index characterization: extreme trend
drought_trend <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/SPEI-6_extreme_trend.tif'))
drought_trend <- -1 * drought_trend # For getting SPEI original units
intensification <- terra::zonal(x = drought_trend, z = wrl, fun = 'mean', na.rm = T); rm(drought_trend)
intensification$Continent <- wrl$Continent
names(intensification)[1] <- 'Average'
gg2 <- intensification |>
  ggplot2::ggplot(aes(x = reorder(Continent, -Average), y = Average)) +
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::xlab('') +
  ggplot2::ylab('SPEI-6 trend 1980-2022') +
  ggplot2::coord_flip() +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'plain', family = 'serif'))

# Drought severity ----
# SPEI-6 time series by growing season
spei6_one_ts <- terra::rast('//CATALOGUE/AgroclimExtremes/agex_indices/agex_spei-6/agex_spei-6_25km/one_s1_spei-6_25km.tif')
spei6_one_ts <- -1 * spei6_one_ts
spei6_two_s1_ts <- terra::rast('//CATALOGUE/AgroclimExtremes/agex_indices/agex_spei-6/agex_spei-6_25km/two_s1_spei-6_25km.tif')
spei6_two_s1_ts <- -1 * spei6_two_s1_ts
spei6_two_s2_ts <- terra::rast('//CATALOGUE/AgroclimExtremes/agex_indices/agex_spei-6/agex_spei-6_25km/two_s2_spei-6_25km.tif')
spei6_two_s2_ts <- -1 * spei6_two_s2_ts

# We choose the more extreme conditions between the growing seasons
spei6_two_ts <- min(spei6_two_s1_ts, spei6_two_s2_ts); rm(spei6_two_s1_ts, spei6_two_s2_ts)
spei6_ts <- terra::merge(spei6_one_ts, spei6_two_ts); rm(spei6_one_ts, spei6_two_ts)

# Time partitions
tp1 <- 1980:1994 # Tertile 1
tp2 <- 1995:2008 # Tertile 2
tp3 <- 2009:2022 # Tertile 3

continents <- c('North America','South America','Africa','Europe','Asia','Oceania')
severity <- continents |>
  purrr::map(.f = function(continent) {
    
    # Time series
    r <- spei6_ts
    r <- r |> terra::crop(wrl[wrl$Continent == continent,]) |> terra::mask(wrl[wrl$Continent == continent,])
    
    # Extreme drought clusters by continent
    agex_sgn_aux <- agex_sgn |> terra::crop(wrl[wrl$Continent == continent,]) |> terra::mask(wrl[wrl$Continent == continent,])
    
    # Number of years with No drought
    tp1_no_drg <- sum(r[[paste0('spei-6_',tp1)]] > -0.5)
    tp2_no_drg <- sum(r[[paste0('spei-6_',tp2)]] > -0.5)
    tp3_no_drg <- sum(r[[paste0('spei-6_',tp3)]] > -0.5)
    tp1_no_smm <- terra::zonal(x = tp1_no_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    tp2_no_smm <- terra::zonal(x = tp2_no_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    tp3_no_smm <- terra::zonal(x = tp3_no_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    names(tp1_no_smm)[2] <- 'NoDrought_1980-1994'
    names(tp2_no_smm)[2] <- 'NoDrought_1995-2008'
    names(tp3_no_smm)[2] <- 'NoDrought_2009-2022'
    no_drought <- Reduce(dplyr::left_join, list(tp1_no_smm, tp2_no_smm, tp3_no_smm))
    # Number of years with Mild drought
    tp1_mi_drg <- sum(r[[paste0('spei-6_',tp1)]] > -1 & r[[paste0('spei-6_',tp1)]] < -0.5)
    tp2_mi_drg <- sum(r[[paste0('spei-6_',tp2)]] > -1 & r[[paste0('spei-6_',tp2)]] < -0.5)
    tp3_mi_drg <- sum(r[[paste0('spei-6_',tp3)]] > -1 & r[[paste0('spei-6_',tp3)]] < -0.5)
    tp1_mi_smm <- terra::zonal(x = tp1_mi_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    tp2_mi_smm <- terra::zonal(x = tp2_mi_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    tp3_mi_smm <- terra::zonal(x = tp3_mi_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    names(tp1_mi_smm)[2] <- 'Mild_1980-1994'
    names(tp2_mi_smm)[2] <- 'Mild_1995-2008'
    names(tp3_mi_smm)[2] <- 'Mild_2009-2022'
    mild <- Reduce(dplyr::left_join, list(tp1_mi_smm, tp2_mi_smm, tp3_mi_smm))
    # Number of years with Moderate drought
    tp1_mo_drg <- sum(r[[paste0('spei-6_',tp1)]] > -1.5 & r[[paste0('spei-6_',tp1)]] < -1)
    tp2_mo_drg <- sum(r[[paste0('spei-6_',tp2)]] > -1.5 & r[[paste0('spei-6_',tp2)]] < -1)
    tp3_mo_drg <- sum(r[[paste0('spei-6_',tp3)]] > -1.5 & r[[paste0('spei-6_',tp3)]] < -1)
    tp1_mo_smm <- terra::zonal(x = tp1_mo_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    tp2_mo_smm <- terra::zonal(x = tp2_mo_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    tp3_mo_smm <- terra::zonal(x = tp3_mo_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    names(tp1_mo_smm)[2] <- 'Moderate_1980-1994'
    names(tp2_mo_smm)[2] <- 'Moderate_1995-2008'
    names(tp3_mo_smm)[2] <- 'Moderate_2009-2022'
    moderate <- Reduce(dplyr::left_join, list(tp1_mo_smm, tp2_mo_smm, tp3_mo_smm))
    # Number of years with Severe drought
    tp1_sv_drg <- sum(r[[paste0('spei-6_',tp1)]] > -2 & r[[paste0('spei-6_',tp1)]] < -1.5)
    tp2_sv_drg <- sum(r[[paste0('spei-6_',tp2)]] > -2 & r[[paste0('spei-6_',tp2)]] < -1.5)
    tp3_sv_drg <- sum(r[[paste0('spei-6_',tp3)]] > -2 & r[[paste0('spei-6_',tp3)]] < -1.5)
    tp1_sv_smm <- terra::zonal(x = tp1_sv_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    tp2_sv_smm <- terra::zonal(x = tp2_sv_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    tp3_sv_smm <- terra::zonal(x = tp3_sv_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    names(tp1_sv_smm)[2] <- 'Severe_1980-1994'
    names(tp2_sv_smm)[2] <- 'Severe_1995-2008'
    names(tp3_sv_smm)[2] <- 'Severe_2009-2022'
    severe <- Reduce(dplyr::left_join, list(tp1_sv_smm, tp2_sv_smm, tp3_sv_smm))
    # Number of years with Extreme drought
    tp1_xt_drg <- sum(r[[paste0('spei-6_',tp1)]] < -2)
    tp2_xt_drg <- sum(r[[paste0('spei-6_',tp2)]] < -2)
    tp3_xt_drg <- sum(r[[paste0('spei-6_',tp3)]] < -2)
    tp1_xt_smm <- terra::zonal(x = tp1_xt_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    tp2_xt_smm <- terra::zonal(x = tp2_xt_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    tp3_xt_smm <- terra::zonal(x = tp3_xt_drg, z = agex_sgn_aux, fun = 'mean', na.rm = T)
    names(tp1_xt_smm)[2] <- 'Extreme_1980-1994'
    names(tp2_xt_smm)[2] <- 'Extreme_1995-2008'
    names(tp3_xt_smm)[2] <- 'Extreme_2009-2022'
    extreme <- Reduce(dplyr::left_join, list(tp1_xt_smm, tp2_xt_smm, tp3_xt_smm))
    
    # Merge counts
    smm <- Reduce(dplyr::left_join, list(extreme, severe, moderate, mild, no_drought))
    dfm <- smm |>
      tidyr::pivot_longer(cols = 2:ncol(smm), names_to = 'Condition', values_to = 'Value') |>
      tidyr::separate(Condition, c('Severity','Time'), sep = '_') |>
      dplyr::mutate(Severity = factor(Severity, levels = c('NoDrought','Mild','Moderate','Severe','Extreme')),
                    Time = factor(Time, levels = c('1980-1994','1995-2008','2009-2022'), ordered = T))
    dfm$Continent <- continent
    return(dfm)
    
}) |> dplyr::bind_rows()

gg3 <- severity |>
  dplyr::mutate(Continent = factor(Continent, levels = continents)) |>
  tidyr::drop_na() |>
  dplyr::group_by(Severity, Time, Continent) |>
  dplyr::summarise(Value = mean(Value)) |>
  dplyr::ungroup() |>
  dplyr::mutate(Severity = factor(Severity, levels = c('NoDrought','Mild','Moderate','Severe','Extreme'), labels = c('No drought','Mild','Moderate','Severe','Extreme'))) |>
  ggplot2::ggplot(aes(x = Time, y = Value, fill = Severity)) +
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::xlab('') +
  ggplot2::ylab('Number of years') +
  ggplot2::scale_fill_manual(values = c('#9dc0cc','#f2d165','#e48451','#a84644','#3e2714')) +
  ggplot2::ylim(c(0, 16)) +
  ggplot2::facet_wrap(~Continent, ncol = 2) +
  ggplot2::coord_flip() +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'bottom',
                 legend.title    = element_blank(),
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 strip.text      = element_text(size = 15, colour = 'black', family = 'serif'),
                 legend.text     = element_text(size = 12, family = 'serif'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'plain', family = 'serif'))

lyt <- matrix(c(1,3,2,3), ncol = 2, nrow = 2, byrow = T)

fg_01 <- ggpubr::annotate_figure(gg1,
                                 top = text_grob(
                                   label  = 'a)',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_02 <- ggpubr::annotate_figure(gg2,
                                 top = text_grob(
                                   label  = 'b)',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_03 <- ggpubr::annotate_figure(gg3,
                                 top = text_grob(
                                   label  = 'c)',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fig2 <- gridExtra::grid.arrange(fg_01, fg_02, fg_03, layout_matrix = lyt)
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/FigureS2-agex-clusters_characterization.png'), plot = fig2, device = 'png', width = 12, height = 10, units = 'in', dpi = 350)
