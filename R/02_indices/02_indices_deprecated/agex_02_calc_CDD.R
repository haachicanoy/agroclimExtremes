## ------------------------------------------ ##
## Maximum number of consecutive dry days
## By: Harold Achicanoy
## WUR & ABC
## Dec. 2023
## ------------------------------------------ ##

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,lubridate))

# Source agro-climatic indices functions
source('https://raw.githubusercontent.com/haachicanoy/agroclimExtremes/main/R/agex_00_indices.R')
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')

# Directories
ddir <- '//CATALOGUE/WFP_ClimateRiskPr1'
root <- '//CATALOGUE/AgroclimExtremes'

# Precipitation files
prc <- list.files(path = paste0(ddir,'/1.Data/AgERA5/precipitation_flux'), pattern = '.nc$', full.names = T)

# ------------------------------------------ #
# One season
# ------------------------------------------ #

# There are two possibilities. One-year growing season which indicates that
# start and ending dates are within the same year. Two-years growing season
# indicating that the start and ending dates cover two consecutive years.

# >>> Processing of One-year growing season
s1_ini <- terra::rast(paste0(root,'/agex_raw_data/agex_one_s1_ini_10km_croplands.tif'))
s1_end <- terra::rast(paste0(root,'/agex_raw_data/agex_one_s1_end_10km_croplands.tif'))
s1_dff <- s1_end - s1_ini
# Selecting the pixels with One-year growing season
s1_dff[s1_end < s1_ini] <- NA
s1_ini[is.na(s1_dff)] <- NA
s1_end[is.na(s1_dff)] <- NA

# >>> Indices computation
# Maximum number of consecutive dry days
lapply(X = yrs, FUN = function(yr){
  prc_flt <- prc[grep(pattern = paste0('_',yr,'[0-9][0-9][0-9][0-9]_'), x = prc)] # Files filtered per year
  rnf <- terra::rast(prc_flt); gc(T)
  rsl <- terra::rapp(x = rnf, first = s1_ini, last = s1_end, fun = cdd); gc(T)
  names(rsl) <- paste0('cdd_',yr)
  terra::writeRaster(x = rsl, filename = paste0(root,'/agex_indices/agex_cdd/agex_cdd_intermediate/one_s1_y1_cdd_',yr,'.tif'), overwrite = T); gc(T)
  return('Done.\n')
})

# >>> Processing of Two-years growing season
s1_ini <- terra::rast(paste0(root,'/agex_raw_data/agex_one_s1_ini_10km_croplands.tif'))
s1_end <- terra::rast(paste0(root,'/agex_raw_data/agex_one_s1_end_10km_croplands.tif'))
s1_dff <- s1_end - s1_ini
# Selecting the pixels with Two-years growing season
s1_dff[s1_end > s1_ini] <- NA
s1_ini[is.na(s1_dff)] <- NA
s1_end[is.na(s1_dff)] <- NA

# >>> Indices computation
# Maximum number of consecutive dry days
lapply(X = yrs, FUN = function(yr){
  # Precipitation files filtered
  if(yr < 2023){
    prc_flt <- prc[unlist(grep2(pattern = c(paste0('_',yr,'[0-9][0-9][0-9][0-9]_'),paste0('_',yr+1,'[0-9][0-9][0-9][0-9]_')), x = prc))] # To implement
    rnf <- terra::rast(prc_flt)
    dys <- ifelse(lubridate::leap_year(yr),366,365)
    rsl <- terra::rapp(x = rnf, first = s1_ini, last = s1_end+dys, fun = cdd)
    terra::writeRaster(x = rsl, filename = paste0(root,'/agex_indices/agex_cdd/agex_cdd_intermediate/one_s1_y2_cdd_',yr,'.tif'), overwrite = T); gc(T)
    return('Done.\n')
  } else {
    break
    return('Done.\n')
  }
})

# >>> Merging results of One-year and Two-years growing seasons
# Maximum number of consecutive dry days
yrs <- 1979:2022
lapply(X = yrs, FUN = function(yr){
  fls <- list.files(pattern = paste0('_',yr,'.tif'), path = paste0(root,'/agex_indices/agex_cdd/agex_cdd_intermediate'), full.names = T)
  cdd <- lapply(fls, terra::rast)
  cdd <- terra::merge(x = cdd[[1]], y = cdd[[2]])
  terra::writeRaster(x = cdd, filename = paste0(root,'/agex_indices/agex_cdd/agex_cdd_10km/one_s1_cdd_',yr,'.tif'), overwrite = T)
})

# ------------------------------------------ #
# Two seasons
# ------------------------------------------ #

# There are two possibilities. One-year growing season which indicates that
# start and ending dates are within the same year. Two-years growing season
# indicating that the start and ending dates cover two consecutive years.
# But in this case, it applies for both growing seasons.

## Season 1
# >>> Processing of One-year growing season
s1_ini <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s1_ini_10km_croplands.tif'))
s1_end <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s1_end_10km_croplands.tif'))
s1_dff <- s1_end - s1_ini
# Selecting the pixels with One-year growing season
s1_dff[s1_end < s1_ini] <- NA
s1_ini[is.na(s1_dff)] <- NA
s1_end[is.na(s1_dff)] <- NA

# >>> Indices computation
# Maximum number of consecutive dry days
lapply(X = yrs, FUN = function(yr){
  prc_flt <- prc[grep(pattern = paste0('_',yr,'[0-9][0-9][0-9][0-9]_'), x = prc)] # Files filtered per year
  rnf <- terra::rast(prc_flt); gc(T)
  rsl <- terra::rapp(x = rnf, first = s1_ini, last = s1_end, fun = cdd); gc(T)
  names(rsl) <- paste0('cdd_',yr)
  terra::writeRaster(x = rsl, filename = paste0(root,'/agex_indices/agex_cdd/agex_cdd_intermediate/two_s1_y1_cdd_',yr,'.tif'), overwrite = T); gc(T)
  return('Done.\n')
})

# >>> Processing of Two-years growing season
s1_ini <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s1_ini_10km_croplands.tif'))
s1_end <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s1_end_10km_croplands.tif'))
s1_dff <- s1_end - s1_ini
# Selecting the pixels with Two-years growing season
s1_dff[s1_end > s1_ini] <- NA
s1_ini[is.na(s1_dff)] <- NA
s1_end[is.na(s1_dff)] <- NA

# >>> Indices computation
# Maximum number of consecutive dry days
lapply(X = yrs, FUN = function(yr){
  # Precipitation files filtered
  if(yr < 2023){
    prc_flt <- prc[unlist(grep2(pattern = c(paste0('_',yr,'[0-9][0-9][0-9][0-9]_'),paste0('_',yr+1,'[0-9][0-9][0-9][0-9]_')), x = prc))] # To implement
    rnf <- terra::rast(prc_flt)
    dys <- ifelse(lubridate::leap_year(yr),366,365)
    rsl <- terra::rapp(x = rnf, first = s1_ini, last = s1_end+dys, fun = cdd)
    names(rsl) <- paste0('cdd_',yr)
    terra::writeRaster(x = rsl, filename = paste0(root,'/agex_indices/agex_cdd/agex_cdd_intermediate/two_s1_y2_cdd_',yr,'.tif'), overwrite = T); gc(T)
    return('Done.\n')
  } else {
    break
    return('Done.\n')
  }
})

# >>> Merging results of One-year and Two-years growing seasons
# Maximum number of consecutive dry days
yrs <- 1979:2022
lapply(X = yrs, FUN = function(yr){
  fls <- list.files(pattern = paste0('two_s1_y[0-9]_cdd_',yr,'.tif'), path = paste0(root,'/agex_indices/agex_cdd/agex_cdd_intermediate'), full.names = T)
  cdd <- lapply(fls, terra::rast)
  cdd <- terra::merge(x = cdd[[1]], y = cdd[[2]])
  terra::writeRaster(x = cdd, filename = paste0(root,'/agex_indices/agex_cdd/agex_cdd_10km/two_s1_cdd_',yr,'.tif'), overwrite = T)
})

## Season 2
# >>> Processing of One-year growing season
s2_ini <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s2_ini_10km_croplands.tif'))
s2_end <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s2_end_10km_croplands.tif'))
s2_dff <- s2_end - s2_ini
# Selecting the pixels with One-year growing season
s2_dff[s2_end < s2_ini] <- NA
s2_ini[is.na(s2_dff)] <- NA
s2_end[is.na(s2_dff)] <- NA

# >>> Indices computation
# Maximum number of consecutive dry days
lapply(X = yrs, FUN = function(yr){
  prc_flt <- prc[grep(pattern = paste0('_',yr,'[0-9][0-9][0-9][0-9]_'), x = prc)] # Files filtered per year
  rnf <- terra::rast(prc_flt); gc(T)
  rsl <- terra::rapp(x = rnf, first = s2_ini, last = s2_end, fun = cdd); gc(T)
  names(rsl) <- paste0('cdd_',yr)
  terra::writeRaster(x = rsl, filename = paste0(root,'/agex_indices/agex_cdd/agex_cdd_intermediate/two_s2_y1_cdd_',yr,'.tif'), overwrite = T); gc(T)
  return('Done.\n')
})

# >>> Processing of Two-years growing season
s2_ini <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s2_ini_10km_croplands.tif'))
s2_end <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s2_end_10km_croplands.tif'))
s2_dff <- s2_end - s2_ini
# Selecting the pixels with Two-years growing season
s2_dff[s2_end > s2_ini] <- NA
s2_ini[is.na(s2_dff)] <- NA
s2_end[is.na(s2_dff)] <- NA

# >>> Indices computation
# Maximum number of consecutive dry days
lapply(X = yrs, FUN = function(yr){
  # Precipitation files filtered
  if(yr < 2023){
    prc_flt <- prc[unlist(grep2(pattern = c(paste0('_',yr,'[0-9][0-9][0-9][0-9]_'),paste0('_',yr+1,'[0-9][0-9][0-9][0-9]_')), x = prc))] # To implement
    rnf <- terra::rast(prc_flt)
    dys <- ifelse(lubridate::leap_year(yr),366,365)
    rsl <- terra::rapp(x = rnf, first = s2_ini, last = s2_end+dys, fun = cdd)
    names(rsl) <- paste0('cdd_',yr)
    terra::writeRaster(x = rsl, filename = paste0(root,'/agex_indices/agex_cdd/agex_cdd_intermediate/two_s2_y2_cdd_',yr,'.tif'), overwrite = T); gc(T)
    return('Done.\n')
  } else {
    break
    return('Done.\n')
  }
})

# >>> Merging results of One-year and Two-years growing seasons
# Maximum number of consecutive dry days
yrs <- 1979:2022
lapply(X = yrs, FUN = function(yr){
  fls <- list.files(pattern = paste0('two_s2_y[0-9]_cdd_',yr,'.tif'), path = paste0(root,'/agex_indices/agex_cdd/agex_cdd_intermediate'), full.names = T)
  cdd <- lapply(fls, terra::rast)
  cdd <- terra::merge(x = cdd[[1]], y = cdd[[2]])
  terra::writeRaster(x = cdd, filename = paste0(root,'/agex_indices/agex_cdd/agex_cdd_10km/two_s2_cdd_',yr,'.tif'), overwrite = T)
})
