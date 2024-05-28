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

### North America ----
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
# One growing season
cls <- 442
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
