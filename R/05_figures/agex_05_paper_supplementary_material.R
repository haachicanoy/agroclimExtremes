## ------------------------------------------ ##
## Supplementary material
## By: Harold Achicanoy
## WUR & ABC
## May 2024
## ------------------------------------------ ##

## R options and packages loading ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,dplyr,stringr,tidyr,pivottabler,data.table,openxlsx))

grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')

## Setup arguments ----
root   <- '//CATALOGUE/WFP_ClimateRiskPr1'
index  <- 'spei-6'

## Extreme drought clusters ----
agex_sgn <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn) <- 'extreme_cluster'

agex_sgn_poly <- terra::vect(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/vct_agex_global_spei-6.gpkg'))

## Arable lands ----
# Cropland areas from MapSPAM 2010
crops_fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_raw_data/croplands/spam2010'), pattern = '_A.tif', full.names = T)

# Crops classification
grp <- read.csv('https://raw.githubusercontent.com/wri/MAPSPAM/master/metadata_tables/4-Methodology-Crops-of-SPAM-2005-2015-02-26.csv')
grp <- grp[,c('SPAM.short.name','GROUP')]
foods <- grp[!(grp$GROUP %in% c('fibres','stimulant','')),'SPAM.short.name']
crops_fls <- crops_fls[grep2(pattern = toupper(foods), x = crops_fls)]

crops_area <- terra::rast(crops_fls) |> terra::mask(mask = agex_sgn_poly)
names(crops_area) <- gsub('spam2010V2r0_global_H_', '', names(crops_area))
names(crops_area) <- gsub('_A', '', names(crops_area))

grp <- read.csv('https://raw.githubusercontent.com/wri/MAPSPAM/master/metadata_tables/4-Methodology-Crops-of-SPAM-2005-2015-02-26.csv')
grp <- grp[,c('SPAM.short.name','SPAM.long.name')]
grp$SPAM.short.name <- toupper(grp$SPAM.short.name)

names(crops_area) <- grp$SPAM.long.name[match(x = names(crops_area), table = grp$SPAM.short.name)]
names(crops_area) <- stringr::str_to_title(names(crops_area))

min_vls <- terra::zonal(x = crops_area, z = agex_sgn_poly, fun = 'min', exact = T, na.rm = T)
min_vls <- cbind(data.frame(extreme_cluster = 1:nrow(min_vls)), min_vls)
max_vls <- terra::zonal(x = crops_area, z = agex_sgn_poly, fun = 'max', exact = T, na.rm = T)
max_vls <- cbind(data.frame(extreme_cluster = 1:nrow(max_vls)), max_vls)
avg_vls <- terra::zonal(x = crops_area, z = agex_sgn_poly, fun = 'mean', exact = T, na.rm = T)
avg_vls <- cbind(data.frame(extreme_cluster = 1:nrow(avg_vls)), avg_vls)
sum_vls <- terra::zonal(x = crops_area, z = agex_sgn_poly, fun = 'sum', na.rm = T)
sum_vls <- cbind(data.frame(extreme_cluster = 1:nrow(sum_vls)), sum_vls)
std_vls <- terra::zonal(x = crops_area, z = agex_sgn_poly, fun = 'sd', na.rm = T)
std_vls <- cbind(data.frame(extreme_cluster = 1:nrow(std_vls)), std_vls)

min_vls <- min_vls |>
  tidyr::pivot_longer(cols = 2:ncol(min_vls), names_to = 'Crop', values_to = 'Min') |>
  base::as.data.frame()
min_vls$Min <- round(min_vls$Min, 2)

max_vls <- max_vls |>
  tidyr::pivot_longer(cols = 2:ncol(max_vls), names_to = 'Crop', values_to = 'Max') |>
  base::as.data.frame()
max_vls$Max <- round(max_vls$Max, 2)

avg_vls <- avg_vls |>
  tidyr::pivot_longer(cols = 2:ncol(avg_vls), names_to = 'Crop', values_to = 'Mean') |>
  base::as.data.frame()
avg_vls$Mean <- round(avg_vls$Mean, 2)

sum_vls <- sum_vls |>
  tidyr::pivot_longer(cols = 2:ncol(sum_vls), names_to = 'Crop', values_to = 'Total') |>
  base::as.data.frame()
sum_vls$Total <- round(sum_vls$Total, 2)

std_vls <- std_vls |>
  tidyr::pivot_longer(cols = 2:ncol(std_vls), names_to = 'Crop', values_to = 'Sd') |>
  base::as.data.frame()
std_vls$Sd <- round(std_vls$Sd, 2)

area_sts_vls <- dplyr::left_join(x = min_vls, y = max_vls, by = c('extreme_cluster','Crop'))
area_sts_vls <- dplyr::left_join(x = area_sts_vls, y = avg_vls, by = c('extreme_cluster','Crop'))
area_sts_vls <- dplyr::left_join(x = area_sts_vls, y = std_vls, by = c('extreme_cluster','Crop'))
area_sts_vls <- dplyr::left_join(x = area_sts_vls, y = sum_vls, by = c('extreme_cluster','Crop'))

# Just totals
totals <- area_sts_vls[,c('extreme_cluster','Crop','Total')] |>
  tidyr::pivot_wider(names_from = 'Crop', values_from = Total) |>
  base::as.data.frame()
write.csv(x = totals, file = paste0(root,'/agroclimExtremes/agex_results/agex_harvested_areas_totals.csv'), row.names = F)

totals_prc <- totals
totals_prc$Total <- rowSums(x = totals[,-1], na.rm = T)
totals_prc[,-1]  <- round(totals_prc[,-1]/totals_prc$Total * 100, 2)
totals_prc$Total <- NULL
write.csv(x = totals_prc, file = paste0(root,'/agroclimExtremes/agex_results/agex_harvested_areas_percentages.csv'), row.names = F)

area_sts_vls_lng <- tidyr::pivot_longer(data = area_sts_vls, cols = `Min`:`Total`, names_to = 'Statistic', values_to = 'Value') |> base::as.data.frame()

pt <- PivotTable$new(evaluationMode = 'batch',
                     processingLibrary = 'data.table')
pt$addData(area_sts_vls_lng)
pt$addColumnDataGroups('Crop', addTotal = F)
pt$addColumnDataGroups('Statistic', addTotal = F)
pt$addRowDataGroups('extreme_cluster', addTotal = F)
pt$defineCalculation(calculationName = 'Summary', type = 'value', valueName = 'Value')
pt$evaluatePivot()

wb <- createWorkbook(creator = Sys.getenv('USERNAME'))
addWorksheet(wb, 'Data')
pt$writeToExcelWorksheet(wb = wb, wsName = 'Data',
                         topRowNumber = 1, leftMostColumnNumber = 1, applyStyles = F)
saveWorkbook(wb, file = paste0(root,'/agroclimExtremes/agex_results/arable_land_statistics_per_crop.xlsx'), overwrite = T)

## Livestock units ----
# Load number of livestock equivalent units
lsu <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_vulnerability/lsu_individuals.tif'))

total_lsu <- terra::zonal(x = lsu, z = agex_sgn_poly, fun = 'sum', na.rm = T)
total_lsu <- round(total_lsu, 2)
total_lsu <- cbind(data.frame(extreme_cluster = 1:nrow(total_lsu)), total_lsu)
write.csv(x = total_lsu, file = paste0(root,'/agroclimExtremes/agex_results/agex_livestock_units_totals.csv'), row.names = F)

prc_lsu <- total_lsu
prc_lsu$Total <- rowSums(prc_lsu[,-1], na.rm = T)
prc_lsu[,-1] <- round(prc_lsu[,-1]/prc_lsu$Total * 100, 2)
prc_lsu$Total <- NULL
write.csv(x = prc_lsu, file = paste0(root,'/agroclimExtremes/agex_results/agex_livestock_units_percentages.csv'), row.names = F)

total_lsu <- total_lsu |>
  tidyr::pivot_longer(cols = 2:ncol(total_lsu), names_to = 'Animal', values_to = 'LSU') |>
  base::as.data.frame()
total_lsu$LSU <- round(total_lsu$LSU, 2)