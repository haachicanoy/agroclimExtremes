## ------------------------------------------ ##
## Compute agro-climatic indices
## By: Harold Achicanoy
## WUR & ABC
## Dec. 2023
## ------------------------------------------ ##

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,lubridate))

# Source agro-climatic indices functions
source('https://raw.githubusercontent.com/haachicanoy/agroclimExtremes/main/agex_0_indices.R')
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')

# Directories
root <- '//CATALOGUE/WFP_ClimateRiskPr1'

# Precipitation files
prc <- list.files(path = paste0(root,'/1.Data/AgERA5/precipitation_flux'), pattern = '.nc$', full.names = T)

# ------------------------------------------ #
# One season
# ------------------------------------------ #
# One-year
s1_ini <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_ini_10km_croplands.tif'))
s1_end <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_end_10km_croplands.tif'))
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
  terra::writeRaster(x = rsl, filename = paste0(root,'/agroclimExtremes/agex_indices/agex_p95/one_s1_y1_p95_',yr,'.tif'), overwrite = T); gc(T)
  return('Done.\n')
})

# Maximum number of consecutive dry days
lapply(X = yrs, FUN = function(yr){
  prc_flt <- prc[grep(pattern = paste0('_',yr,'[0-9][0-9][0-9][0-9]_'), x = prc)] # Files filtered per year
  rnf <- terra::rast(prc_flt); gc(T)
  rsl <- terra::rapp(x = rnf, first = s1_ini, last = s1_end, fun = cdd); gc(T)
  names(rsl) <- paste0('cdd_',yr)
  terra::writeRaster(x = rsl, filename = paste0(root,'/agroclimExtremes/agex_indices/agex_cdd/one_s1_y1_cdd_',yr,'.tif'), overwrite = T); gc(T)
  return('Done.\n')
})

# Two-years
s1_ini <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_ini_10km_croplands.tif'))
s1_end <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_end_10km_croplands.tif'))
s1_dff <- s1_end - s1_ini
s1_dff[s1_end > s1_ini] <- NA
s1_ini[is.na(s1_dff)] <- NA
s1_end[is.na(s1_dff)] <- NA

# Percentile 95th of daily precipitation
yrs <- 1979:2023
lapply(X = yrs, FUN = function(yr){
  # Precipitation files filtered
  if(yr < 2023){
    prc_flt <- prc[unlist(grep2(pattern = c(paste0('_',yr,'[0-9][0-9][0-9][0-9]_'),paste0('_',yr+1,'[0-9][0-9][0-9][0-9]_')), x = prc))] # To implement
    rnf <- terra::rast(prc_flt)
    dys <- ifelse(lubridate::leap_year(yr),366,365)
    rsl <- terra::rapp(x = rnf, first = s1_ini, last = s1_end+dys, fun = p95)
    terra::writeRaster(x = rsl, filename = paste0(root,'/agroclimExtremes/agex_indices/agex_p95/one_s1_y2_p95_',yr,'.tif'), overwrite = T); gc(T)
    return('Done.\n')
  } else {
    break
    return('Done.\n')
  }
})

# Maximum number of consecutive dry days
lapply(X = yrs, FUN = function(yr){
  # Precipitation files filtered
  if(yr < 2023){
    prc_flt <- prc[unlist(grep2(pattern = c(paste0('_',yr,'[0-9][0-9][0-9][0-9]_'),paste0('_',yr+1,'[0-9][0-9][0-9][0-9]_')), x = prc))] # To implement
    rnf <- terra::rast(prc_flt)
    dys <- ifelse(lubridate::leap_year(yr),366,365)
    rsl <- terra::rapp(x = rnf, first = s1_ini, last = s1_end+dys, fun = cdd)
    terra::writeRaster(x = rsl, filename = paste0(root,'/agroclimExtremes/agex_indices/agex_cdd/one_s1_y2_cdd_',yr,'.tif'), overwrite = T); gc(T)
    return('Done.\n')
  } else {
    break
    return('Done.\n')
  }
})

# Merge results
yrs <- 1979:2022
lapply(X = yrs, FUN = function(yr){
  fls <- list.files(pattern = paste0('_',yr,'.tif'), path = paste0(root,'/agroclimExtremes/agex_indices/agex_cdd'), full.names = T)
  cdd <- lapply(fls, terra::rast)
  cdd <- terra::merge(x = cdd[[1]], y = cdd[[2]])
  terra::writeRaster(x = cdd, filename = paste0(root,'/agroclimExtremes/agex_indices/agex_cdd/one_s1_cdd_',yr,'.tif'), overwrite = T)
})

# ------------------------------------------ #
# Two seasons
# ------------------------------------------ #
