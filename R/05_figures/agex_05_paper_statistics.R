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

## Clusters presence ----
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
# Extreme drought clusters in more than one country (%)
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

## Top countries with more than clusters due to their size ----
sort(table(agex_sgn_metrics$countries[agex_sgn_metrics$countries_count == 1]), decreasing = T)[1:5]

## Transcontinental clusters ----
# Absolute frequency
trscntn_cnt <- strsplit(x = agex_sgn_metrics$continents, split = ',') |>
  purrr::map(length) |>
  base::unlist() |>
  table()
# Relative frequency
trscntn_prc <- trscntn_cnt/nrow(agex_sgn_metrics) * 100
round(sum(trscntn_prc[-1]), 1)

## Continental clusters distribution ----
trscntn_id <- strsplit(x = agex_sgn_metrics$continents, split = ',') |>
  purrr::map(length) |>
  base::unlist()
round(sort(table(agex_sgn_metrics$continents[trscntn_id == 1]), decreasing = T)/nrow(agex_sgn_metrics) * 100, 1)

## Potential cooperation per continent ----
agex_sgn_metrics |>
  dplyr::filter(trscntn_id == 1) |>
  dplyr::group_by(continents) |>
  dplyr::summarize(number_of_clusters = dplyr::n(),
                   collaboration = sum(countries_count > 1)) |>
  dplyr::ungroup() |>
  dplyr::mutate(freq = round(collaboration/number_of_clusters * 100, 1)) |>
  dplyr::arrange(-freq) |>
  base::as.data.frame()

## Cooperation case studies per continent ----
table(agex_sgn_metrics$countries_count[grep(pattern = '[aA]frica', agex_sgn_metrics$continents)])/101 * 100
table(agex_sgn_metrics$countries_count[grep(pattern = '[eE]urope', agex_sgn_metrics$continents)])/76 * 100
table(agex_sgn_metrics$countries_count[grep(pattern = '[sS]outh [aA]merica', agex_sgn_metrics$continents)])/77 * 100
table(agex_sgn_metrics$countries_count[grep(pattern = '[nN]orth [aA]merica', agex_sgn_metrics$continents)])/64 * 100
table(agex_sgn_metrics$countries_count[grep(pattern = '[aA]sia', agex_sgn_metrics$continents)])/158 * 100
table(agex_sgn_metrics$countries_count[grep(pattern = '[oO]ceania', agex_sgn_metrics$continents)])/15
