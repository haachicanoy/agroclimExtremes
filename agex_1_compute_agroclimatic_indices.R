## ------------------------------------------ ##
## Compute agro-climatic indices
## By: Harold Achicanoy
## WUR & ABC
## Dec. 2023
## ------------------------------------------ ##

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra))

# Source agro-climatic indices functions
source('https://raw.githubusercontent.com/haachicanoy/agroclimExtremes/main/agex_0_indices.R')

# Directories
root <- '//CATALOGUE/WFP_ClimateRiskPr1'

# Precipitation files
prc <- list.files(path = paste0(root,'/1.Data/AgERA5/precipitation_flux'), pattern = '.nc$', full.names = T)

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
yrs <- 1979:2023
lapply(X = yrs, FUN = function(yr){
  prc_flt <- prc[grep(pattern = paste0('_',yr,'[0-9][0-9][0-9][0-9]_'), x = prc)] # Files filtered per year
  rnf <- terra::rast(prc_flt); gc(T)
  rsl <- terra::rapp(x = rnf, first = s1_ini, last = s1_end, fun = p95); gc(T)
  names(rsl) <- c('p95','day')
  terra::writeRaster(x = rsl, filename = paste0(root,'/agroclimExtremes/indices/p95/one_s1_p95_',yr,'.tif'), overwrite = T); gc(T)
  return('Done.\n')
})

# Maximum number of consecutive dry days
lapply(X = yrs, FUN = function(yr){
  prc_flt <- prc[grep(pattern = paste0('_',yr,'[0-9][0-9][0-9][0-9]_'), x = prc)] # Files filtered per year
  rnf <- terra::rast(prc_flt); gc(T)
  rsl <- terra::rapp(x = rnf, first = s1_ini, last = s1_end, fun = cdd); gc(T)
  names(rsl) <- paste0('cdd_',yr)
  terra::writeRaster(x = rsl, filename = paste0(root,'/agroclimExtremes/indices/cdd/one_s1_cdd_',yr,'.tif'), overwrite = T); gc(T)
  return('Done.\n')
})

p95_idx <- index |> purrr::map(1) |> terra::rast()
p95_whn <- index |> purrr::map(2) |> terra::rast()

# Two-years (TO DO)
yrs <- 1981:2023
p95_index <- lapply(X = yrs, FUN = function(yr){
  # Precipitation files filtered
  prc_flt <- prc[grep(pattern = c(yr,yr+1), x = prc)] # To implement
  rnf <- terra::rast(prc_flt)
  rsl <- terra::rapp(x = rnf, first = s1_ini_5km, last = s1_end_5km, fun = p95)
  terra::writeRaster(x = rsl, filename = paste0(root,'/yr_',yr,'_world.tif'))
  return('Done.\n')
})

# ------------------------------------------ #
# Two seasons
# ------------------------------------------ #
