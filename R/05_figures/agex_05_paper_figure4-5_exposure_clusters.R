# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Figure 4: Crop composition within drought hotspots
# Table S4: Median crop distribution within diversity hotspots by continent
# Figure 5: Livestock unit distribution within drought hotspots
# Table S5: Median livestock unit distribution within economic hotspots by continent
# By: Harold Achicanoy
# WUR & ABC
# Created in 2025
# Modified in February 2026
# --------------------------------------------------------------- #

# R options and user-defined functions ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,corrplot,MetBrewer,sf,rnaturalearth,tidyverse,ggpubr,tseries,biscale,pals,cowplot,gridExtra,grid))

capitalize_first_word <- function(text_string) {
  # Check if the input is a character string
  if (!is.character(text_string) || length(text_string) == 0) {
    return(text_string) # Return as is if not a valid string
  }
  
  # Capitalize the first letter
  first_letter <- toupper(substr(text_string, 1, 1))
  
  # Get the rest of the string
  rest_of_string <- substr(text_string, 2, nchar(text_string))
  
  # Combine them
  paste0(first_letter, rest_of_string)
}

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

# Crops and livestock diversity ----
if (!exists('intensification_vs_cropdiversity')) {
  # Crop classes diversity
  crp_diversity <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/crop_classes_diversity_25km.tif'))
  
  # Index characterization: extreme trend
  idx_sxt <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/SPEI-6_extreme_trend.tif'))
  
  # Zonal statistics per extreme drought cluster (mean)
  agex_sgn_extreme_trend <- terra::zonal(x = idx_sxt, z = agex_sgn_num, fun = 'mean', na.rm = T, as.raster = T)
  metrics1 <- terra::zonal(x = idx_sxt, z = agex_sgn_num, fun = 'mean', na.rm = T)
  agex_sgn_crp_diversity <- terra::zonal(x = crp_diversity, z = agex_sgn_num, fun = 'mean', na.rm = T, as.raster = T)
  metrics2 <- terra::zonal(x = crp_diversity, z = agex_sgn_num, fun = 'mean', na.rm = T)
  metrics <- dplyr::left_join(x = metrics1, y = metrics2, by = 'extreme_cluster'); rm(metrics1, metrics2)
  
  # Merge data.frames
  agex_sgn_severity_vs_crops <- terra::as.data.frame(x = c(agex_sgn, agex_sgn_extreme_trend, agex_sgn_crp_diversity), xy = T, cell = T)
  rm(agex_sgn_extreme_trend, agex_sgn_crp_diversity)
  
  # Produce bivariate maps
  # Define severity and diversity quantiles
  severity_brks <- quantile(x = metrics$`SPEI-6_slope_95th`, probs = seq(0,1,1/4))
  severity_brks[1] <- terra::as.data.frame(idx_sxt) |> min()
  severity_brks[5] <- terra::as.data.frame(idx_sxt) |> max()
  diversity_brks <- quantile(x = metrics$crop_classes_diversity, probs = seq(0,1,1/4))
  diversity_brks[1] <- terra::as.data.frame(crp_diversity)|> min()
  diversity_brks[5] <- terra::as.data.frame(crp_diversity)|> max()
  # Create categories
  agex_sgn_severity_vs_crops$Severity_class <- cut(x = agex_sgn_severity_vs_crops$`SPEI-6_slope_95th`, breaks = severity_brks) |> as.numeric() |> as.character()
  agex_sgn_severity_vs_crops$Diversity_class <- cut(x = agex_sgn_severity_vs_crops$crop_classes_diversity, breaks = diversity_brks) |> as.numeric() |> as.character()
  # Create bivariate categories
  agex_sgn_severity_vs_crops$bi_class <- paste0(agex_sgn_severity_vs_crops$Severity_class,'-',agex_sgn_severity_vs_crops$Diversity_class)
  agex_sgn_severity_vs_crops <- agex_sgn_severity_vs_crops[,c('cell','x','y','extreme_cluster','SPEI-6_slope_95th','crop_classes_diversity','bi_class')]
  # Define color palette for bivariate categories
  colours <- data.frame(bi_class = sort(unique(agex_sgn_severity_vs_crops$bi_class)))
  colours$category <- 1:nrow(colours)
  agex_sgn_severity_vs_crops <- dplyr::left_join(x = agex_sgn_severity_vs_crops, y = colours, by = 'bi_class')
  
  intensification_vs_cropdiversity <- agex_sgn_severity_vs_crops[,c('extreme_cluster','bi_class')] |> unique()
  rm(colours, agex_sgn_severity_vs_crops, diversity_brks, severity_brks, metrics, idx_sxt, crp_diversity)
}
names(intensification_vs_cropdiversity)[2] <- 'Drought_vs_CropDiversity'
intensification_vs_cropdiversity$extreme_cluster <- as.integer(intensification_vs_cropdiversity$extreme_cluster)

if (!exists('intensification_vs_lvstdiversity')) {
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
  
  intensification_vs_lvstdiversity <- agex_sgn_severity_vs_lvst[,c('extreme_cluster','bi_class')] |> unique()
  rm(colours, agex_sgn_severity_vs_lvst, diversity_brks, severity_brks, metrics, idx_sxt, lsu_diversity)
}
names(intensification_vs_lvstdiversity)[2] <- 'Drought_vs_LivestockDiversity'
intensification_vs_lvstdiversity$extreme_cluster <- as.integer(intensification_vs_lvstdiversity$extreme_cluster)

# Crops and livestock economic value ----
if (!exists('intensification_vs_cropvop')) {
  # Crop VOP
  crp_vop <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/food_vop.tif'))
  
  # Index characterization: extreme trend
  idx_sxt <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/SPEI-6_extreme_trend.tif'))
  
  # Zonal statistics per extreme drought cluster (mean)
  agex_sgn_extreme_trend <- terra::zonal(x = idx_sxt, z = agex_sgn_num, fun = 'mean', na.rm = T, as.raster = T)
  metrics1 <- terra::zonal(x = idx_sxt, z = agex_sgn_num, fun = 'mean', na.rm = T)
  agex_sgn_crp_vop <- terra::zonal(x = crp_vop, z = agex_sgn_num, fun = 'mean', na.rm = T, as.raster = T)
  metrics2 <- terra::zonal(x = crp_vop, z = agex_sgn_num, fun = 'mean', na.rm = T)
  metrics <- dplyr::left_join(x = metrics1, y = metrics2, by = 'extreme_cluster'); rm(metrics1, metrics2)
  
  # Merge data.frames
  agex_sgn_severity_vs_crops <- terra::as.data.frame(x = c(agex_sgn, agex_sgn_extreme_trend, agex_sgn_crp_vop), xy = T, cell = T)
  rm(agex_sgn_extreme_trend, agex_sgn_crp_vop)
  
  # Produce bivariate maps
  # Define severity and diversity quantiles
  severity_brks <- quantile(x = metrics$`SPEI-6_slope_95th`, probs = seq(0,1,1/4))
  severity_brks[1] <- terra::as.data.frame(idx_sxt) |> min()
  severity_brks[5] <- terra::as.data.frame(idx_sxt) |> max()
  diversity_brks <- quantile(x = metrics$food_vop, probs = seq(0,1,1/4))
  diversity_brks[1] <- terra::as.data.frame(crp_vop)|> min()
  diversity_brks[5] <- terra::as.data.frame(crp_vop)|> max()
  # Create categories
  agex_sgn_severity_vs_crops$Severity_class <- cut(x = agex_sgn_severity_vs_crops$`SPEI-6_slope_95th`, breaks = severity_brks) |> as.numeric() |> as.character()
  agex_sgn_severity_vs_crops$Diversity_class <- cut(x = agex_sgn_severity_vs_crops$food_vop, breaks = diversity_brks) |> as.numeric() |> as.character()
  # Create bivariate categories
  agex_sgn_severity_vs_crops$bi_class <- paste0(agex_sgn_severity_vs_crops$Severity_class,'-',agex_sgn_severity_vs_crops$Diversity_class)
  agex_sgn_severity_vs_crops <- agex_sgn_severity_vs_crops[,c('cell','x','y','extreme_cluster','SPEI-6_slope_95th','food_vop','bi_class')]
  # Define color palette for bivariate categories
  colours <- data.frame(bi_class = sort(unique(agex_sgn_severity_vs_crops$bi_class)))
  colours$category <- 1:nrow(colours)
  agex_sgn_severity_vs_crops <- dplyr::left_join(x = agex_sgn_severity_vs_crops, y = colours, by = 'bi_class')
  
  intensification_vs_cropvop <- agex_sgn_severity_vs_crops[,c('extreme_cluster','bi_class')] |> unique()
  rm(colours, agex_sgn_severity_vs_crops, diversity_brks, severity_brks, metrics, idx_sxt, crp_vop)
}
names(intensification_vs_cropvop)[2] <- 'Drought_vs_CropVoP'
intensification_vs_cropvop$extreme_cluster <- as.integer(intensification_vs_cropvop$extreme_cluster)

if (!exists('intensification_vs_lvstunits')) {
  # Livestock equivalent units
  lvstck_units  <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/lsu_total_25km.tif'))
  
  # Index characterization: extreme trend
  idx_sxt <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/SPEI-6_extreme_trend.tif'))
  
  # Zonal statistics per extreme drought cluster (mean)
  agex_sgn_extreme_trend <- terra::zonal(x = idx_sxt, z = agex_sgn_num, fun = 'mean', na.rm = T, as.raster = T)
  metrics1 <- terra::zonal(x = idx_sxt, z = agex_sgn_num, fun = 'mean', na.rm = T)
  agex_sgn_lvstck_units <- terra::zonal(x = lvstck_units, z = agex_sgn_num, fun = 'mean', na.rm = T, as.raster = T)
  metrics2 <- terra::zonal(x = lvstck_units, z = agex_sgn_num, fun = 'mean', na.rm = T)
  metrics <- dplyr::left_join(x = metrics1, y = metrics2, by = 'extreme_cluster'); rm(metrics1, metrics2)
  
  # Merge data.frames
  agex_sgn_severity_vs_lvstck <- terra::as.data.frame(x = c(agex_sgn, agex_sgn_extreme_trend, agex_sgn_lvstck_units), xy = T, cell = T)
  rm(agex_sgn_extreme_trend, agex_sgn_lvstck_units)
  
  # Produce bivariate maps
  # Define severity and diversity quantiles
  severity_brks <- quantile(x = metrics$`SPEI-6_slope_95th`, probs = seq(0,1,1/4))
  severity_brks[1] <- terra::as.data.frame(idx_sxt) |> min()
  severity_brks[5] <- terra::as.data.frame(idx_sxt) |> max()
  diversity_brks <- quantile(x = metrics$total_livestock_units, probs = seq(0,1,1/4))
  diversity_brks[1] <- terra::as.data.frame(lvstck_units)|> min()
  diversity_brks[5] <- terra::as.data.frame(lvstck_units)|> max()
  # Create categories
  agex_sgn_severity_vs_lvstck$Severity_class <- cut(x = agex_sgn_severity_vs_lvstck$`SPEI-6_slope_95th`, breaks = severity_brks) |> as.numeric() |> as.character()
  agex_sgn_severity_vs_lvstck$Diversity_class <- cut(x = agex_sgn_severity_vs_lvstck$total_livestock_units, breaks = diversity_brks) |> as.numeric() |> as.character()
  # Create bivariate categories
  agex_sgn_severity_vs_lvstck$bi_class <- paste0(agex_sgn_severity_vs_lvstck$Severity_class,'-',agex_sgn_severity_vs_lvstck$Diversity_class)
  agex_sgn_severity_vs_lvstck <- agex_sgn_severity_vs_lvstck[,c('cell','x','y','extreme_cluster','SPEI-6_slope_95th','total_livestock_units','bi_class')]
  # Define color palette for bivariate categories
  colours <- data.frame(bi_class = sort(unique(agex_sgn_severity_vs_lvstck$bi_class)))
  colours$category <- 1:nrow(colours)
  agex_sgn_severity_vs_lvstck <- dplyr::left_join(x = agex_sgn_severity_vs_lvstck, y = colours, by = 'bi_class')
  
  intensification_vs_lvstunits <- agex_sgn_severity_vs_lvstck[,c('extreme_cluster','bi_class')] |> unique()
  rm(colours, agex_sgn_severity_vs_lvstck, diversity_brks, severity_brks, metrics, idx_sxt, lvstck_units)
}
names(intensification_vs_lvstunits)[2] <- 'Drought_vs_LivestockUnits'
intensification_vs_lvstunits$extreme_cluster <- as.integer(intensification_vs_lvstunits$extreme_cluster)

hotspot_categories <- c('4-4')

intensification_vs_cropdiversity |> dplyr::filter(Drought_vs_CropDiversity %in% hotspot_categories)
intensification_vs_lvstdiversity |> dplyr::filter(Drought_vs_LivestockDiversity %in% hotspot_categories)
intensification_vs_cropvop |> dplyr::filter(Drought_vs_CropVoP %in% hotspot_categories)
intensification_vs_lvstunits |> dplyr::filter(Drought_vs_LivestockUnits %in% hotspot_categories)

agex_sgn_metrics <- dplyr::left_join(x = agex_sgn_metrics, y = intensification_vs_cropdiversity, by = 'extreme_cluster')
agex_sgn_metrics <- dplyr::left_join(x = agex_sgn_metrics, y = intensification_vs_cropvop, by = 'extreme_cluster')
agex_sgn_metrics <- dplyr::left_join(x = agex_sgn_metrics, y = intensification_vs_lvstdiversity, by = 'extreme_cluster')
agex_sgn_metrics <- dplyr::left_join(x = agex_sgn_metrics, y = intensification_vs_lvstunits, by = 'extreme_cluster')

rm(intensification_vs_cropdiversity, intensification_vs_lvstdiversity, intensification_vs_cropvop, intensification_vs_lvstunits)

View(agex_sgn_metrics)

# Hotspots counts and overlap ----
# Food security hotspots
table((agex_sgn_metrics$Drought_vs_CropDiversity == '4-4') +
        (agex_sgn_metrics$Drought_vs_LivestockDiversity == '4-4'))

# Economic value hotspots
table((agex_sgn_metrics$Drought_vs_CropVoP == '4-4') +
        (agex_sgn_metrics$Drought_vs_LivestockUnits == '4-4'))

grep(pattern = '4-|3-', x = agex_sgn_metrics$Drought_vs_CropDiversity) |>length()
nrow(agex_sgn_metrics)

verify <- grep(pattern = 'Asia', x = agex_sgn_metrics$continents)
length(verify)
continent_subdata <- agex_sgn_metrics[verify,]
continent_subdata[grep(pattern = '4-|3-', x = continent_subdata$Drought_vs_CropDiversity),] |> dim()

sum(continent_subdata$Drought_vs_LivestockDiversity == '4-4')
continent_subdata[continent_subdata$Drought_vs_LivestockDiversity == '4-4',]

table((continent_subdata$Drought_vs_CropDiversity == '4-4') +
        (continent_subdata$Drought_vs_LivestockDiversity == '4-4'))
table((continent_subdata$Drought_vs_CropVoP == '4-4') +
        (continent_subdata$Drought_vs_LivestockUnits == '4-4'))

# Crop counts
crop_fls <- list.files(path = file.path(root,'agex_raw_data'), pattern = '*count_10km.tif$', full.names = T)
crop_nms <- basename(crop_fls) |> gsub(pattern = 'agex_cropclass_', replacement = '', x = _) |> gsub(pattern = '_count_10km.tif', replacement = '', x = _)
crop_cnt <- terra::rast(crop_fls); rm(crop_fls)
names(crop_cnt) <- crop_nms; rm(crop_nms)

# Livestock units
lsu <- terra::rast(file.path(root,'agex_results/agex_vulnerability/lsu_individuals.tif'))

Threatened_CropDiversity <- agex_sgn_metrics |>
  dplyr::filter(Drought_vs_CropDiversity == '4-4') |>
  dplyr::pull(extreme_cluster) |>
  purrr::map(.f = function(xtr_clus) {
    
    # Extreme cluster of interest
    aux <- agex_sgn
    aux <- (aux == xtr_clus) * 1
    aux <- terra::classify(x = aux, cbind(0, NA))
    crd <- terra::as.data.frame(aux, xy = T, na.rm = T)
    xtd <- terra::ext(c(min(crd$x), max(crd$x), min(crd$y), max(crd$y))); rm(crd)
    aux <- terra::crop(x = aux, y = xtd)
    vct <- terra::as.polygons(aux); rm(aux,crd,xtd)
    
    crop_classes_count <- crop_cnt
    crop_classes_count <- crop_classes_count |> terra::crop(terra::ext(vct)) |> terra::mask(vct)
    
    crops_count <- summary(crop_classes_count) |>
      as.data.frame() |>
      tidyr::separate(Freq, into = c('Statistics','Value'), sep = ':') |>
      dplyr::select(Var2, Statistics, Value) |>
      dplyr::filter(Statistics == 'Median ') |>
      dplyr::mutate(Value = as.numeric(Value),
                    Value = round(Value, 0))
    names(crops_count)[1] <- 'Crop'
    crops_count$extreme_cluster <- xtr_clus
    
    return(crops_count)
    
  }) |> dplyr::bind_rows()

Threatened_CropDiversity$Statistics <- NULL
tcd <- Threatened_CropDiversity |> tidyr::pivot_wider(names_from = 'Crop', values_from = Value) |> base::as.data.frame()
tcd <- dplyr::left_join(x = agex_sgn_metrics[agex_sgn_metrics$Drought_vs_CropDiversity == '4-4',c('extreme_cluster','growing_seasons','countries_count','isos','countries','continents')], y = tcd, by = 'extreme_cluster')

# --------------------------------------------------------------- #
# Figure 4
# --------------------------------------------------------------- #
gg1 <- dplyr::left_join(x = Threatened_CropDiversity,
                 y = agex_sgn_metrics[agex_sgn_metrics$Drought_vs_CropDiversity == '4-4',c('extreme_cluster','growing_seasons','countries_count','isos','countries','continents')],
                 by = 'extreme_cluster') |>
  dplyr::group_by(Crop, continents) |>
  dplyr::summarise(Value = median(Value)) |>
  dplyr::ungroup() |>
  base::as.data.frame() |>
  dplyr::filter(continents != 'Asia,Europe') |>
  dplyr::mutate(Crop = trimws(Crop, which = 'left'),
                Crop = gsub(pattern = '\\_|\\.', replacement = ' ', x = Crop),
                Crop = capitalize_first_word(Crop)) |>
  ggplot2::ggplot(aes(x = reorder(continents, Value), y = Value, fill = Crop)) +
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::xlab('') +
  ggplot2::ylab('Median number of crop by class') +
  ggplot2::coord_flip() +
  ggplot2::theme_minimal() +
  ggplot2::scale_fill_manual(values = c('#332288','#88CCEE','#44AA99','#CC6677','#999933',
                                        '#DDDDDD','#DDCC77','#882255','#AA4499','#117733')) + # MetBrewer::met.brewer(name = 'Signac', n = 10)
  ggplot2::theme(legend.position = 'bottom',
                 legend.title    = element_blank(),
                 axis.text       = element_text(size = 19, family = 'serif'),
                 axis.text.x     = element_text(size = 19, colour = 'black'),
                 axis.text.y     = element_text(size = 19, colour = 'black'),
                 axis.title.x    = element_text(size = 19, colour = 'black'),
                 legend.text     = element_text(size = 15, family = 'serif'),
                 axis.title      = element_text(size = 19, colour = 'black', face = 'plain', family = 'serif'))
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/Figure4_paper1-crop_diversity_exposure.png'), plot = gg1, device = 'png', width = 11, height = 8, units = 'in', dpi = 350)

# Table S4
dplyr::left_join(x = Threatened_CropDiversity,
                 y = agex_sgn_metrics[agex_sgn_metrics$Drought_vs_CropDiversity == '4-4',c('extreme_cluster','growing_seasons','countries_count','isos','countries','continents')],
                 by = 'extreme_cluster') |>
  dplyr::group_by(Crop, continents) |>
  dplyr::summarise(Value = median(Value)) |>
  dplyr::ungroup() |>
  base::as.data.frame() |>
  dplyr::filter(continents != 'Asia,Europe') |>
  dplyr::mutate(Crop = trimws(Crop, which = 'left'),
                Crop = gsub(pattern = '\\_|\\.', replacement = ' ', x = Crop),
                Crop = capitalize_first_word(Crop)) |>
  tidyr::pivot_wider(names_from = 'continents', values_from = 'Value') |>
  base::as.data.frame() |>
  write.csv(x = _, file = 'D:/median_number_crops_in_hotspots.csv', row.names = F)

# --------------------------------------------------------------- #
# --------------------------------------------------------------- #

Threatened_LivestockUnits <- agex_sgn_metrics |>
  dplyr::filter(Drought_vs_LivestockUnits == '4-4') |>
  dplyr::pull(extreme_cluster) |>
  purrr::map(.f = function(xtr_clus) {
    
    # Extreme cluster of interest
    aux <- agex_sgn
    aux <- (aux == xtr_clus) * 1
    aux <- terra::classify(x = aux, cbind(0, NA))
    crd <- terra::as.data.frame(aux, xy = T, na.rm = T)
    xtd <- terra::ext(c(min(crd$x), max(crd$x), min(crd$y), max(crd$y))); rm(crd)
    aux <- terra::crop(x = aux, y = xtd)
    vct <- terra::as.polygons(aux); rm(aux,crd,xtd)
    
    livestock_units <- lsu
    livestock_units <- livestock_units |> terra::crop(terra::ext(vct)) |> terra::mask(vct)
    
    animals_count <- summary(livestock_units) |>
      as.data.frame() |>
      tidyr::separate(Freq, into = c('Statistics','Value'), sep = ':') |>
      dplyr::select(Var2, Statistics, Value) |>
      dplyr::filter(Statistics == 'Median ') |>
      dplyr::mutate(Value = as.numeric(Value),
                    Value = round(Value, 0))
    names(animals_count)[1] <- 'Animal'
    animals_count$extreme_cluster <- xtr_clus
    
    return(animals_count)
    
  }) |> dplyr::bind_rows()

Threatened_LivestockUnits$Statistics <- NULL
tlu <- Threatened_LivestockUnits |> tidyr::pivot_wider(names_from = 'Animal', values_from = Value) |> base::as.data.frame()
tlu <- dplyr::left_join(x = agex_sgn_metrics[agex_sgn_metrics$Drought_vs_LivestockUnits == '4-4',c('extreme_cluster','growing_seasons','countries_count','isos','countries','continents')], y = tlu, by = 'extreme_cluster')

# --------------------------------------------------------------- #
# Figure 5
# --------------------------------------------------------------- #
gg2 <- dplyr::left_join(x = Threatened_LivestockUnits,
                 y = agex_sgn_metrics[agex_sgn_metrics$Drought_vs_LivestockUnits == '4-4',c('extreme_cluster','growing_seasons','countries_count','isos','countries','continents')],
                 by = 'extreme_cluster') |>
  dplyr::group_by(Animal, continents) |>
  dplyr::summarise(Value = median(Value)) |>
  dplyr::ungroup() |>
  base::as.data.frame() |>
  dplyr::mutate(Animal = trimws(Animal, which = 'left')) |>
  ggplot2::ggplot(aes(x = reorder(continents, Value), y = Value, fill = Animal)) +
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::xlab('') +
  ggplot2::ylab('Livestock units median') +
  ggplot2::coord_flip() +
  ggplot2::theme_minimal() +
  ggplot2::scale_fill_manual(values = c('#332288','#88CCEE','#44AA99','#117733',
                                        '#DDCC77','#CC6677','#882255')) + # MetBrewer::met.brewer(name = 'Signac', n = 7)
  ggplot2::theme(legend.position = 'bottom',
                 legend.title    = element_blank(),
                 axis.text       = element_text(size = 19, family = 'serif'),
                 axis.text.x     = element_text(size = 19, colour = 'black'),
                 axis.text.y     = element_text(size = 19, colour = 'black'),
                 axis.title.x    = element_text(size = 19, colour = 'black'),
                 legend.text     = element_text(size = 15, family = 'serif'),
                 axis.title      = element_text(size = 19, colour = 'black', face = 'plain', family = 'serif'))
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/Figure4_paper1-livestock_units_exposure.png'), plot = gg2, device = 'png', width = 11, height = 8, units = 'in', dpi = 350)

# --------------------------------------------------------------- #
# Table S5
# --------------------------------------------------------------- #
dplyr::left_join(x = Threatened_LivestockUnits,
                 y = agex_sgn_metrics[agex_sgn_metrics$Drought_vs_LivestockUnits == '4-4',c('extreme_cluster','growing_seasons','countries_count','isos','countries','continents')],
                 by = 'extreme_cluster') |>
  dplyr::group_by(Animal, continents) |>
  dplyr::summarise(Value = median(Value)) |>
  dplyr::ungroup() |>
  base::as.data.frame() |>
  dplyr::mutate(Animal = trimws(Animal, which = 'left')) |>
  tidyr::pivot_wider(names_from = 'continents', values_from = 'Value') |>
  base::as.data.frame() |>
  write.csv(x = _, file = 'D:/median_number_LSU_in_hotspots.csv', row.names = F)

# agex_sgn_metrics |>
#   dplyr::select(continents, area, total_population) |>
#   dplyr::group_by(continents) |>
#   dplyr::summarise(nclusters = n(),
#                    area_avg    = mean(area),
#                    pop_avg     = mean(total_population),
#                    area_sd     = sd(area, na.rm = T),
#                    pop_sd      = sd(total_population, na.rm = T)) |>
#   dplyr::filter(continents %in% c('North America','South America','Africa','Europe','Asia','Oceania')) |>
#   base::as.data.frame() |>
#   View()
# 
# agex_sgn_metrics |>
#   dplyr::filter(Drought_vs_CropVoP == '4-4' | Drought_vs_LivestockUnits == '4-4') |>
#   write.csv(x = _, file = 'D:/OneDrive - CGIAR/PhD/papers/paper1/submission/Achicanoy_et_al-SM-Exposure_hotspots.csv', row.names = F)
