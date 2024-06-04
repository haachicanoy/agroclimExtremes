## ------------------------------------------ ##
## Paper statistics
## By: Harold Achicanoy
## WUR & ABC
## May 2024
## ------------------------------------------ ##

## R options and packages loading ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,MetBrewer,rnaturalearth,tidyverse,ggpubr,tseries,biscale,pals,cowplot,gridExtra,grid))

# Activate Windows fonts
extrafont::font_import()
extrafont::loadfonts(device = 'win')
extrafont::fonts()

## Setup arguments ----
root   <- '//CATALOGUE/WFP_ClimateRiskPr1'
index  <- 'spei-6'

## Extreme drought clusters ----
# Categorical
agex_sgn <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn) <- 'extreme_cluster'
cls <- data.frame(extreme_cluster = sort(unique(terra::values(agex_sgn))))
cls$id <- cls$extreme_cluster
cls <- cls[,c('id','extreme_cluster')]
cls$extreme_cluster <- as.character(cls$extreme_cluster)
levels(agex_sgn) <- cls; rm(cls)

# Numerical
agex_sgn_num <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn_num) <- 'extreme_cluster'

# Metrics
agex_sgn_metrics <- utils::read.csv(paste0(root,'/agroclimExtremes/agex_results/agex_all_metrics.csv'))

## Extreme drought clusters ----

### Clusters presence ----
cat('Extreme drought clusters in one country (%):',
    round(table(agex_sgn_metrics$countries_count)[[1]]/sum(table(agex_sgn_metrics$countries_count)) * 100, 1),
    '\n')
cat('Extreme drought clusters in two or more countries (%):',
    100 - round(table(agex_sgn_metrics$countries_count)[[1]]/sum(table(agex_sgn_metrics$countries_count)) * 100, 1),
    '\n')

## One growing season
cat('One-GS. Extreme drought clusters in one country (%):',
    round(table(agex_sgn_metrics$countries_count[agex_sgn_metrics$growing_seasons == 1])[[1]]/sum(table(agex_sgn_metrics$countries_count[agex_sgn_metrics$growing_seasons == 1])) * 100, 1),
    '\n')
cat('One-GS. Extreme drought clusters in two or more countries (%):',
    100 - round(table(agex_sgn_metrics$countries_count[agex_sgn_metrics$growing_seasons == 1])[[1]]/sum(table(agex_sgn_metrics$countries_count[agex_sgn_metrics$growing_seasons == 1])) * 100, 1),
    '\n')

## Two growing seasons
cat('Two-GS. Extreme drought clusters in one country (%):',
    round(table(agex_sgn_metrics$countries_count[agex_sgn_metrics$growing_seasons == 2])[[1]]/sum(table(agex_sgn_metrics$countries_count[agex_sgn_metrics$growing_seasons == 2])) * 100, 1),
    '\n')
cat('Two-GS. Extreme drought clusters in two or more countries (%):',
    100 - round(table(agex_sgn_metrics$countries_count[agex_sgn_metrics$growing_seasons == 2])[[1]]/sum(table(agex_sgn_metrics$countries_count[agex_sgn_metrics$growing_seasons == 2])) * 100, 2),
    '\n')

### Top countries with more than clusters due to their size ----
sort(table(agex_sgn_metrics$countries[agex_sgn_metrics$countries_count == 1]), decreasing = T)[1:5]

### Transcontinental clusters ----
# Absolute frequency
trscntn_cnt <- strsplit(x = agex_sgn_metrics$continents, split = ',') |>
  purrr::map(length) |>
  base::unlist() |>
  table()
# Relative frequency
trscntn_prc <- trscntn_cnt/nrow(agex_sgn_metrics) * 100
round(sum(trscntn_prc[-1]), 1)

### Continental clusters distribution ----
trscntn_id <- strsplit(x = agex_sgn_metrics$continents, split = ',') |>
  purrr::map(length) |>
  base::unlist()
round(sort(table(agex_sgn_metrics$continents[trscntn_id == 1]), decreasing = T)/nrow(agex_sgn_metrics) * 100, 1)

### Potential cooperation per continent ----
agex_sgn_metrics |>
  dplyr::filter(trscntn_id == 1) |>
  dplyr::group_by(continents) |>
  dplyr::summarize(number_of_clusters = dplyr::n(),
                   collaboration = sum(countries_count > 1)) |>
  dplyr::ungroup() |>
  dplyr::mutate(freq = round(collaboration/number_of_clusters * 100, 1)) |>
  dplyr::arrange(-freq) |>
  base::as.data.frame()

### Cooperation case studies per continent ----
table(agex_sgn_metrics$countries_count[grep(pattern = '[aA]frica', agex_sgn_metrics$continents)])/101 * 100
table(agex_sgn_metrics$countries_count[grep(pattern = '[eE]urope', agex_sgn_metrics$continents)])/76 * 100
table(agex_sgn_metrics$countries_count[grep(pattern = '[sS]outh [aA]merica', agex_sgn_metrics$continents)])/77 * 100
table(agex_sgn_metrics$countries_count[grep(pattern = '[nN]orth [aA]merica', agex_sgn_metrics$continents)])/64 * 100
table(agex_sgn_metrics$countries_count[grep(pattern = '[aA]sia', agex_sgn_metrics$continents)])/158 * 100
table(agex_sgn_metrics$countries_count[grep(pattern = '[oO]ceania', agex_sgn_metrics$continents)])/15

## Drought extremes and agricultural exposure hotspots ----
# Identification of vulnerable areas
severity_brks <- quantile(x = agex_sgn_metrics$`SPEI.6_slope_95th`, probs = seq(0,1,1/4))
diversity_brks <- quantile(x = agex_sgn_metrics$crop_classes_diversity, probs = seq(0,1,1/4))
# Create categories
agex_sgn_metrics$Severity_class <- cut(x = agex_sgn_metrics$`SPEI.6_slope_95th`, breaks = severity_brks) |> as.numeric() |> as.character()
agex_sgn_metrics$Diversity_class <- cut(x = agex_sgn_metrics$crop_classes_diversity, breaks = diversity_brks) |> as.numeric() |> as.character()

agex_sgn_metrics$agricultural_exposure <- paste0(agex_sgn_metrics$Severity_class,'-',agex_sgn_metrics$Diversity_class)
agex_sgn_metrics$Severity_class <- NULL
agex_sgn_metrics$Diversity_class <- NULL

severity_brks <- quantile(x = agex_sgn_metrics$`SPEI.6_slope_95th`, probs = seq(0,1,1/4))
diversity_brks <- quantile(x = agex_sgn_metrics$livestock_units_total, probs = seq(0,1,1/4))

agex_sgn_metrics$Severity_class <- cut(x = agex_sgn_metrics$`SPEI.6_slope_95th`, breaks = severity_brks) |> as.numeric() |> as.character()
agex_sgn_metrics$Diversity_class <- cut(x = agex_sgn_metrics$livestock_units_total, breaks = diversity_brks) |> as.numeric() |> as.character()

agex_sgn_metrics$livestock_exposure <- paste0(agex_sgn_metrics$Severity_class,'-',agex_sgn_metrics$Diversity_class)
agex_sgn_metrics$Severity_class <- NULL
agex_sgn_metrics$Diversity_class <- NULL
rm(severity_brks, diversity_brks)

mtx_ref <- expand.grid(1:4, 1:4)
names(mtx_ref) <- c('Severity','Diversity')
mtx_ref <- mtx_ref |>
  dplyr::arrange(Severity) |>
  dplyr::mutate(key = paste0(Severity,'-',Diversity)) |>
  base::as.data.frame()

### North America ----
NA_crops <- agex_sgn_metrics$agricultural_exposure[grep(pattern = 'North America', x = agex_sgn_metrics$continents)] |>
  gtools::mixedsort() |>
  table() |>
  base::as.data.frame() |>
  tidyr::separate(col = 'Var1', into = c('Severity','Diversity'), remove = F)
names(NA_crops)[1] <- 'key'
NA_crops$key <- as.character(NA_crops$key)
mss_cats <- data.frame(key = base::setdiff(x = mtx_ref$key, y = NA_crops$key))
if(nrow(mss_cats) > 0){
  mss_cats <- mss_cats |>
    tidyr::separate(col = 'key', into = c('Severity','Diversity'), remove = F) |>
    dplyr::mutate(Freq = 0) |>
    base::as.data.frame()
  NA_crops <- rbind(NA_crops, mss_cats); rm(mss_cats)
  NA_crops <- NA_crops |>
    dplyr::arrange(Severity, Diversity) |>
    base::as.data.frame()
}
NA_crops |>
  dplyr::filter(Severity >= 3 & Diversity >= 3)

NA_lvstc <- agex_sgn_metrics$livestock_exposure[grep(pattern = 'North America', x = agex_sgn_metrics$continents)] |>
  gtools::mixedsort() |>
  table() |>
  base::as.data.frame() |>
  tidyr::separate(col = 'Var1', into = c('Severity','Diversity'), remove = F)
names(NA_lvstc)[1] <- 'key'
NA_lvstc$key <- as.character(NA_lvstc$key)
mss_cats <- data.frame(key = base::setdiff(x = mtx_ref$key, y = NA_lvstc$key))
if(nrow(mss_cats) > 0){
  mss_cats <- mss_cats |>
    tidyr::separate(col = 'key', into = c('Severity','Diversity'), remove = F) |>
    dplyr::mutate(Freq = 0) |>
    base::as.data.frame()
  NA_lvstc <- rbind(NA_lvstc, mss_cats); rm(mss_cats)
  NA_lvstc <- NA_lvstc |>
    dplyr::arrange(Severity, Diversity) |>
    base::as.data.frame()
}
NA_lvstc |>
  dplyr::filter(Severity >= 3 & Diversity >= 3)

# One growing season
cls <- 223
# What crops are grown in a specific extreme cluster
arable_land_sts <- read.csv(paste0(root,'/agroclimExtremes/agex_results/agex_harvested_areas_percentages.csv'))
arable_land_sts |>
  dplyr::filter(extreme_cluster == cls) |>
  tidyr::pivot_longer(-1, names_to = 'Crop', values_to = 'Percentage') |>
  dplyr::arrange(-Percentage)

# During when (one month before, both cases)
s1_ini <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_ini_25km_croplands.tif'))
names(s1_ini) <- 'S1_start'
s1_end <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_end_25km_croplands.tif'))
names(s1_end) <- 'S1_end'

clndr <- terra::zonal(x = c(s1_ini, s1_end), z = agex_sgn, fun = 'median', na.rm = T) |>
  dplyr::filter(extreme_cluster == cls) |>
  base::as.data.frame()

clndr <- rbind(clndr,
               data.frame(extreme_cluster = cls,
                          S1_start = lubridate::month(as.Date(clndr$S1_start,'2014-01-01')),
                          S1_end = lubridate::month(as.Date(clndr$S1_end,'2014-01-01'))
               )
)
clndr$Date_format <- c('DOY','Month')

# How many people
agex_sgn_metrics[agex_sgn_metrics$extreme_cluster == cls, 'total_population']

## How much agricultural value vs total
# Value of production (food crops)
agex_sgn_metrics$food_vop[agex_sgn_metrics$extreme_cluster == cls]
# Value of production (total)
agex_sgn_metrics$total_vop[agex_sgn_metrics$extreme_cluster == cls]
# Percentage
round(agex_sgn_metrics$food_vop[agex_sgn_metrics$extreme_cluster == cls]/agex_sgn_metrics$total_vop[agex_sgn_metrics$extreme_cluster == cls] * 100, 1)

### South America ----
SA_crops <- agex_sgn_metrics$agricultural_exposure[grep(pattern = 'South America', x = agex_sgn_metrics$continents)] |>
  gtools::mixedsort() |>
  table() |>
  base::as.data.frame() |>
  tidyr::separate(col = 'Var1', into = c('Severity','Diversity'), remove = F)
names(SA_crops)[1] <- 'key'
SA_crops$key <- as.character(SA_crops$key)
mss_cats <- data.frame(key = base::setdiff(x = mtx_ref$key, y = SA_crops$key))
if(nrow(mss_cats) > 0){
  mss_cats <- mss_cats |>
    tidyr::separate(col = 'key', into = c('Severity','Diversity'), remove = F) |>
    dplyr::mutate(Freq = 0) |>
    base::as.data.frame()
  SA_crops <- rbind(SA_crops, mss_cats); rm(mss_cats)
  SA_crops <- SA_crops |>
    dplyr::arrange(Severity, Diversity) |>
    base::as.data.frame()
}
SA_crops |>
  dplyr::filter(Severity >= 3 & Diversity >= 3)

SA_lvstc <- agex_sgn_metrics$livestock_exposure[grep(pattern = 'South America', x = agex_sgn_metrics$continents)] |>
  gtools::mixedsort() |>
  table() |>
  base::as.data.frame() |>
  tidyr::separate(col = 'Var1', into = c('Severity','Diversity'), remove = F)
names(SA_lvstc)[1] <- 'key'
SA_lvstc$key <- as.character(SA_lvstc$key)
mss_cats <- data.frame(key = base::setdiff(x = mtx_ref$key, y = SA_lvstc$key))
if(nrow(mss_cats) > 0){
  mss_cats <- mss_cats |>
    tidyr::separate(col = 'key', into = c('Severity','Diversity'), remove = F) |>
    dplyr::mutate(Freq = 0) |>
    base::as.data.frame()
  SA_lvstc <- rbind(SA_lvstc, mss_cats); rm(mss_cats)
  SA_lvstc <- SA_lvstc |>
    dplyr::arrange(Severity, Diversity) |>
    base::as.data.frame()
}
SA_lvstc |>
  dplyr::filter(Severity >= 3 & Diversity >= 3)

# Two growing seasons
cls <- 420
# What crops are grown in a specific extreme cluster
arable_land_sts <- read.csv(paste0(root,'/agroclimExtremes/agex_results/agex_harvested_areas_percentages.csv'))
arable_land_sts |>
  dplyr::filter(extreme_cluster == cls) |>
  tidyr::pivot_longer(-1, names_to = 'Crop', values_to = 'Percentage') |>
  dplyr::arrange(-Percentage)

# During when (one month before, both cases)
s1_ini <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_two_s1_ini_25km_croplands.tif'))
names(s1_ini) <- 'S1_start'
s1_end <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_two_s1_end_25km_croplands.tif'))
names(s1_end) <- 'S1_end'
s2_ini <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_two_s2_ini_25km_croplands.tif'))
names(s2_ini) <- 'S2_start'
s2_end <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_two_s2_end_25km_croplands.tif'))
names(s2_end) <- 'S2_end'

clndr <- terra::zonal(x = c(s1_ini, s1_end, s2_ini, s2_end), z = agex_sgn, fun = 'median', na.rm = T) |>
  dplyr::filter(extreme_cluster == cls) |>
  base::as.data.frame()
clndr <- rbind(clndr,
               data.frame(extreme_cluster = cls,
                          S1_start = lubridate::month(as.Date(clndr$S1_start,'2014-01-01')),
                          S1_end = lubridate::month(as.Date(clndr$S1_end,'2014-01-01')),
                          S2_start = lubridate::month(as.Date(clndr$S2_start,'2014-01-01')),
                          S2_end = lubridate::month(as.Date(clndr$S2_end,'2014-01-01'))
                          )
               )
clndr$Date_format <- c('DOY','Month')

# How many people
agex_sgn_metrics[agex_sgn_metrics$extreme_cluster == cls, 'total_population']

## How much agricultural value vs total
# Value of production (food crops)
agex_sgn_metrics$food_vop[agex_sgn_metrics$extreme_cluster == cls]
# Value of production (total)
agex_sgn_metrics$total_vop[agex_sgn_metrics$extreme_cluster == cls]
# Percentage
round(agex_sgn_metrics$food_vop[agex_sgn_metrics$extreme_cluster == cls]/agex_sgn_metrics$total_vop[agex_sgn_metrics$extreme_cluster == cls] * 100, 1)

### Africa ----
AF_crops <- agex_sgn_metrics$agricultural_exposure[grep(pattern = 'Africa', x = agex_sgn_metrics$continents)] |>
  gtools::mixedsort() |>
  table() |>
  base::as.data.frame() |>
  tidyr::separate(col = 'Var1', into = c('Severity','Diversity'), remove = F)
names(AF_crops)[1] <- 'key'
AF_crops$key <- as.character(AF_crops$key)
mss_cats <- data.frame(key = base::setdiff(x = mtx_ref$key, y = AF_crops$key))
if(nrow(mss_cats) > 0){
  mss_cats <- mss_cats |>
    tidyr::separate(col = 'key', into = c('Severity','Diversity'), remove = F) |>
    dplyr::mutate(Freq = 0) |>
    base::as.data.frame()
  AF_crops <- rbind(AF_crops, mss_cats); rm(mss_cats)
  AF_crops <- AF_crops |>
    dplyr::arrange(Severity, Diversity) |>
    base::as.data.frame()
}
AF_crops |>
  dplyr::filter(Severity >= 3 & Diversity >= 3)

AF_lvstc <- agex_sgn_metrics$livestock_exposure[grep(pattern = 'Africa', x = agex_sgn_metrics$continents)] |>
  gtools::mixedsort() |>
  table() |>
  base::as.data.frame() |>
  tidyr::separate(col = 'Var1', into = c('Severity','Diversity'), remove = F)
names(AF_lvstc)[1] <- 'key'
AF_lvstc$key <- as.character(AF_lvstc$key)
mss_cats <- data.frame(key = base::setdiff(x = mtx_ref$key, y = AF_lvstc$key))
if(nrow(mss_cats) > 0){
  mss_cats <- mss_cats |>
    tidyr::separate(col = 'key', into = c('Severity','Diversity'), remove = F) |>
    dplyr::mutate(Freq = 0) |>
    base::as.data.frame()
  AF_lvstc <- rbind(AF_lvstc, mss_cats); rm(mss_cats)
  AF_lvstc <- AF_lvstc |>
    dplyr::arrange(Severity, Diversity) |>
    base::as.data.frame()
}
AF_lvstc |>
  dplyr::filter(Severity >= 3 & Diversity >= 3)

# One growing season
cls <- 270
# What crops are grown in a specific extreme cluster
arable_land_sts <- read.csv(paste0(root,'/agroclimExtremes/agex_results/agex_harvested_areas_percentages.csv'))
arable_land_sts |>
  dplyr::filter(extreme_cluster == cls) |>
  tidyr::pivot_longer(-1, names_to = 'Crop', values_to = 'Percentage') |>
  dplyr::arrange(-Percentage)

# During when (one month before, both cases)
s1_ini <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_ini_25km_croplands.tif'))
names(s1_ini) <- 'S1_start'
s1_end <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_end_25km_croplands.tif'))
names(s1_end) <- 'S1_end'

clndr <- terra::zonal(x = c(s1_ini, s1_end), z = agex_sgn, fun = 'median', na.rm = T) |>
  dplyr::filter(extreme_cluster == cls) |>
  base::as.data.frame()

clndr <- rbind(clndr,
               data.frame(extreme_cluster = cls,
                          S1_start = lubridate::month(as.Date(clndr$S1_start,'2014-01-01')),
                          S1_end = lubridate::month(as.Date(clndr$S1_end,'2014-01-01'))
               )
)
clndr$Date_format <- c('DOY','Month')

# How many people
agex_sgn_metrics[agex_sgn_metrics$extreme_cluster == cls, 'total_population']

## How much agricultural value vs total
# Value of production (food crops)
agex_sgn_metrics$food_vop[agex_sgn_metrics$extreme_cluster == cls]
# Value of production (total)
agex_sgn_metrics$total_vop[agex_sgn_metrics$extreme_cluster == cls]
# Percentage
round(agex_sgn_metrics$food_vop[agex_sgn_metrics$extreme_cluster == cls]/agex_sgn_metrics$total_vop[agex_sgn_metrics$extreme_cluster == cls] * 100, 1)

### Europe ----
EU_crops <- agex_sgn_metrics$agricultural_exposure[grep(pattern = 'Europe', x = agex_sgn_metrics$continents)] |>
  gtools::mixedsort() |>
  table() |>
  base::as.data.frame() |>
  tidyr::separate(col = 'Var1', into = c('Severity','Diversity'), remove = F)
names(EU_crops)[1] <- 'key'
EU_crops$key <- as.character(EU_crops$key)
mss_cats <- data.frame(key = base::setdiff(x = mtx_ref$key, y = EU_crops$key))
if(nrow(mss_cats) > 0){
  mss_cats <- mss_cats |>
    tidyr::separate(col = 'key', into = c('Severity','Diversity'), remove = F) |>
    dplyr::mutate(Freq = 0) |>
    base::as.data.frame()
  EU_crops <- rbind(EU_crops, mss_cats); rm(mss_cats)
  EU_crops <- EU_crops |>
    dplyr::arrange(Severity, Diversity) |>
    base::as.data.frame()
}
EU_crops |>
  dplyr::filter(Severity >= 3 & Diversity >= 3)

EU_lvstc <- agex_sgn_metrics$livestock_exposure[grep(pattern = 'Europe', x = agex_sgn_metrics$continents)] |>
  gtools::mixedsort() |>
  table() |>
  base::as.data.frame() |>
  tidyr::separate(col = 'Var1', into = c('Severity','Diversity'), remove = F)
names(EU_lvstc)[1] <- 'key'
EU_lvstc$key <- as.character(EU_lvstc$key)
mss_cats <- data.frame(key = base::setdiff(x = mtx_ref$key, y = EU_lvstc$key))
if(nrow(mss_cats) > 0){
  mss_cats <- mss_cats |>
    tidyr::separate(col = 'key', into = c('Severity','Diversity'), remove = F) |>
    dplyr::mutate(Freq = 0) |>
    base::as.data.frame()
  EU_lvstc <- rbind(EU_lvstc, mss_cats); rm(mss_cats)
  EU_lvstc <- EU_lvstc |>
    dplyr::arrange(Severity, Diversity) |>
    base::as.data.frame()
}
EU_lvstc |>
  dplyr::filter(Severity >= 3 & Diversity >= 3)

# One growing season
cls <- 33
# What crops are grown in a specific extreme cluster
arable_land_sts <- read.csv(paste0(root,'/agroclimExtremes/agex_results/agex_harvested_areas_percentages.csv'))
arable_land_sts |>
  dplyr::filter(extreme_cluster == cls) |>
  tidyr::pivot_longer(-1, names_to = 'Crop', values_to = 'Percentage') |>
  dplyr::arrange(-Percentage)

# During when (one month before, both cases)
s1_ini <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_ini_25km_croplands.tif'))
names(s1_ini) <- 'S1_start'
s1_end <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_end_25km_croplands.tif'))
names(s1_end) <- 'S1_end'

clndr <- terra::zonal(x = c(s1_ini, s1_end), z = agex_sgn, fun = 'median', na.rm = T) |>
  dplyr::filter(extreme_cluster == cls) |>
  base::as.data.frame()

clndr <- rbind(clndr,
               data.frame(extreme_cluster = cls,
                          S1_start = lubridate::month(as.Date(clndr$S1_start,'2014-01-01')),
                          S1_end = lubridate::month(as.Date(clndr$S1_end,'2014-01-01'))
               )
)
clndr$Date_format <- c('DOY','Month')

# How many people
agex_sgn_metrics[agex_sgn_metrics$extreme_cluster == cls, 'total_population']

## How much agricultural value vs total
# Value of production (food crops)
agex_sgn_metrics$food_vop[agex_sgn_metrics$extreme_cluster == cls]
# Value of production (total)
agex_sgn_metrics$total_vop[agex_sgn_metrics$extreme_cluster == cls]
# Percentage
round(agex_sgn_metrics$food_vop[agex_sgn_metrics$extreme_cluster == cls]/agex_sgn_metrics$total_vop[agex_sgn_metrics$extreme_cluster == cls] * 100, 1)

### Asia ----
AS_crops <- agex_sgn_metrics$agricultural_exposure[grep(pattern = 'Asia', x = agex_sgn_metrics$continents)] |>
  gtools::mixedsort() |>
  table() |>
  base::as.data.frame() |>
  tidyr::separate(col = 'Var1', into = c('Severity','Diversity'), remove = F)
names(AS_crops)[1] <- 'key'
AS_crops$key <- as.character(AS_crops$key)
mss_cats <- data.frame(key = base::setdiff(x = mtx_ref$key, y = AS_crops$key))
if(nrow(mss_cats) > 0){
  mss_cats <- mss_cats |>
    tidyr::separate(col = 'key', into = c('Severity','Diversity'), remove = F) |>
    dplyr::mutate(Freq = 0) |>
    base::as.data.frame()
  AS_crops <- rbind(AS_crops, mss_cats); rm(mss_cats)
  AS_crops <- AS_crops |>
    dplyr::arrange(Severity, Diversity) |>
    base::as.data.frame()
}
AS_crops |>
  dplyr::filter(Severity >= 3 & Diversity >= 3)

AS_lvstc <- agex_sgn_metrics$livestock_exposure[grep(pattern = 'Asia', x = agex_sgn_metrics$continents)] |>
  gtools::mixedsort() |>
  table() |>
  base::as.data.frame() |>
  tidyr::separate(col = 'Var1', into = c('Severity','Diversity'), remove = F)
names(AS_lvstc)[1] <- 'key'
AS_lvstc$key <- as.character(AS_lvstc$key)
mss_cats <- data.frame(key = base::setdiff(x = mtx_ref$key, y = AS_lvstc$key))
if(nrow(mss_cats) > 0){
  mss_cats <- mss_cats |>
    tidyr::separate(col = 'key', into = c('Severity','Diversity'), remove = F) |>
    dplyr::mutate(Freq = 0) |>
    base::as.data.frame()
  AS_lvstc <- rbind(AS_lvstc, mss_cats); rm(mss_cats)
  AS_lvstc <- AS_lvstc |>
    dplyr::arrange(Severity, Diversity) |>
    base::as.data.frame()
}
AS_lvstc |>
  dplyr::filter(Severity >= 3 & Diversity >= 3)

# One growing season
cls <- 42
# What crops are grown in a specific extreme cluster
arable_land_sts <- read.csv(paste0(root,'/agroclimExtremes/agex_results/agex_harvested_areas_percentages.csv'))
arable_land_sts |>
  dplyr::filter(extreme_cluster == cls) |>
  tidyr::pivot_longer(-1, names_to = 'Crop', values_to = 'Percentage') |>
  dplyr::arrange(-Percentage)

# During when (one month before, both cases)
s1_ini <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_ini_25km_croplands.tif'))
names(s1_ini) <- 'S1_start'
s1_end <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_end_25km_croplands.tif'))
names(s1_end) <- 'S1_end'

clndr <- terra::zonal(x = c(s1_ini, s1_end), z = agex_sgn, fun = 'median', na.rm = T) |>
  dplyr::filter(extreme_cluster == cls) |>
  base::as.data.frame()

clndr <- rbind(clndr,
               data.frame(extreme_cluster = cls,
                          S1_start = lubridate::month(as.Date(clndr$S1_start,'2014-01-01')),
                          S1_end = lubridate::month(as.Date(clndr$S1_end,'2014-01-01'))
               )
)
clndr$Date_format <- c('DOY','Month')

# How many people
agex_sgn_metrics[agex_sgn_metrics$extreme_cluster == cls, 'total_population']

## How much agricultural value vs total
# Value of production (food crops)
agex_sgn_metrics$food_vop[agex_sgn_metrics$extreme_cluster == cls]
# Value of production (total)
agex_sgn_metrics$total_vop[agex_sgn_metrics$extreme_cluster == cls]
# Percentage
round(agex_sgn_metrics$food_vop[agex_sgn_metrics$extreme_cluster == cls]/agex_sgn_metrics$total_vop[agex_sgn_metrics$extreme_cluster == cls] * 100, 1)

### Oceania ----
OC_crops <- agex_sgn_metrics$agricultural_exposure[grep(pattern = 'Oceania', x = agex_sgn_metrics$continents)] |>
  gtools::mixedsort() |>
  table() |>
  base::as.data.frame() |>
  tidyr::separate(col = 'Var1', into = c('Severity','Diversity'), remove = F)
names(OC_crops)[1] <- 'key'
OC_crops$key <- as.character(OC_crops$key)
mss_cats <- data.frame(key = base::setdiff(x = mtx_ref$key, y = OC_crops$key))
if(nrow(mss_cats) > 0){
  mss_cats <- mss_cats |>
    tidyr::separate(col = 'key', into = c('Severity','Diversity'), remove = F) |>
    dplyr::mutate(Freq = 0) |>
    base::as.data.frame()
  OC_crops <- rbind(OC_crops, mss_cats); rm(mss_cats)
  OC_crops <- OC_crops |>
    dplyr::arrange(Severity, Diversity) |>
    base::as.data.frame()
}
OC_crops |>
  dplyr::filter(Severity >= 2 & Diversity >= 2)

OC_lvstc <- agex_sgn_metrics$livestock_exposure[grep(pattern = 'Oceania', x = agex_sgn_metrics$continents)] |>
  gtools::mixedsort() |>
  table() |>
  base::as.data.frame() |>
  tidyr::separate(col = 'Var1', into = c('Severity','Diversity'), remove = F)
names(OC_lvstc)[1] <- 'key'
OC_lvstc$key <- as.character(OC_lvstc$key)
mss_cats <- data.frame(key = base::setdiff(x = mtx_ref$key, y = OC_lvstc$key))
if(nrow(mss_cats) > 0){
  mss_cats <- mss_cats |>
    tidyr::separate(col = 'key', into = c('Severity','Diversity'), remove = F) |>
    dplyr::mutate(Freq = 0) |>
    base::as.data.frame()
  OC_lvstc <- rbind(OC_lvstc, mss_cats); rm(mss_cats)
  OC_lvstc <- OC_lvstc |>
    dplyr::arrange(Severity, Diversity) |>
    base::as.data.frame()
}
OC_lvstc |>
  dplyr::filter(Severity >= 3 & Diversity >= 3)

# One growing season
cls <- 444
# What crops are grown in a specific extreme cluster
arable_land_sts <- read.csv(paste0(root,'/agroclimExtremes/agex_results/agex_harvested_areas_percentages.csv'))
arable_land_sts |>
  dplyr::filter(extreme_cluster == cls) |>
  tidyr::pivot_longer(-1, names_to = 'Crop', values_to = 'Percentage') |>
  dplyr::arrange(-Percentage)

# During when (one month before, both cases)
s1_ini <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_ini_25km_croplands.tif'))
names(s1_ini) <- 'S1_start'
s1_end <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_end_25km_croplands.tif'))
names(s1_end) <- 'S1_end'

clndr <- terra::zonal(x = c(s1_ini, s1_end), z = agex_sgn, fun = 'median', na.rm = T) |>
  dplyr::filter(extreme_cluster == cls) |>
  base::as.data.frame()

clndr <- rbind(clndr,
               data.frame(extreme_cluster = cls,
                          S1_start = lubridate::month(as.Date(clndr$S1_start,'2014-01-01')),
                          S1_end = lubridate::month(as.Date(clndr$S1_end,'2014-01-01'))
               )
)
clndr$Date_format <- c('DOY','Month')

# How many people
agex_sgn_metrics[agex_sgn_metrics$extreme_cluster == cls, 'total_population']

## How much agricultural value vs total
# Value of production (food crops)
agex_sgn_metrics$food_vop[agex_sgn_metrics$extreme_cluster == cls]
# Value of production (total)
agex_sgn_metrics$total_vop[agex_sgn_metrics$extreme_cluster == cls]
# Percentage
round(agex_sgn_metrics$food_vop[agex_sgn_metrics$extreme_cluster == cls]/agex_sgn_metrics$total_vop[agex_sgn_metrics$extreme_cluster == cls] * 100, 1)
