# ------------------------------------------ #
# Compute agro-climatic indices
# By: Harold Achicanoy
# WUR & ABC
# Dec. 2023
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra))

# Source agro-climatic indices functions
source('https://raw.githubusercontent.com/haachicanoy/agroclimExtremes/main/agex_0_indices.R')

# Directories
root <- '//CATALOGUE/WFP_ClimateRiskPr1'

# Precipitation files
prc <- list.files(path = '//CATALOGUE/Workspace14/WFP_ClimateRiskPr/1.Data/ERA5/precipitation_flux', pattern = '.nc$', full.names = T)

# ------------------------------------------ #
# One season
# ------------------------------------------ #
# One-year
s1_ini <- terra::rast(paste0(root,'/agroclimExtremes/data/results/one_s1_ini_10km.tif'))
s1_end <- terra::rast(paste0(root,'/agroclimExtremes/data/results/one_s1_end_10km.tif'))
s1_dff <- s1_end - s1_ini
s1_dff[s1_end < s1_ini] <- NA
s1_ini[is.na(s1_dff)] <- NA
s1_end[is.na(s1_dff)] <- NA

# Percentile 95th of daily precipitation
yrs <- 1981:2019
lapply(X = yrs, FUN = function(yr){
  prc_flt <- prc[grep(pattern = yr, x = prc)] # Files filtered per year
  rnf <- terra::rast(prc_flt); gc(T)
  rsl <- terra::rapp(x = rnf, first = s1_ini, last = s1_end, fun = p95); gc(T)
  terra::writeRaster(x = rsl, filename = paste0(root,'/one_s1_p95_',yr,'.tif')); gc(T)
  return('Done.\n')
})

# Maximum number of consecutive dry days
lapply(X = yrs, FUN = function(yr){
  prc_flt <- prc[grep(pattern = yr, x = prc)] # Files filtered per year
  rnf <- terra::rast(prc_flt); gc(T)
  rsl <- terra::rapp(x = rnf, first = s1_ini, last = s1_end, fun = cdd); gc(T)
  terra::writeRaster(x = rsl, filename = paste0(root,'/one_s1_cdd_',yr,'.tif')); gc(T)
  return('Done.\n')
})

p95_idx <- index |> purrr::map(1) |> terra::rast()
p95_whn <- index |> purrr::map(2) |> terra::rast()

# Two-years
yrs <- 1981:2023
p95_index <- lapply(X = yrs, FUN = function(yr){
  # Precipitation files filtered
  prc_flt <- prc[grep(pattern = c(yr,yr+1), x = prc)] # To implement
  rnf <- terra::rast(prc_flt)
  rsl <- terra::rapp(x = rnf, first = s1_ini_5km, last = s1_end_5km, fun = p95)
  terra::writeRaster(x = rsl, filename = paste0(root,'/yr_',yr,'_world.tif'), gdal=c('COMPRESS=NONE', 'TFW=YES'), datatype='INT1U')
  return('Done.\n')
})

# ------------------------------------------ #
# Two seasons
# ------------------------------------------ #
