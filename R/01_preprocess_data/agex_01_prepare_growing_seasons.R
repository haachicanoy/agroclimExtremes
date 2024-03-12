# ------------------------------------------ #
# Prepare growing seasons - croplands
# By: Harold Achicanoy
# WUR & ABC
# Dec. 2023
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
rm(list = ls()); gc(TRUE)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra))

# Template rasters
tmp_10km <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5.tif')
tmp_25km <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5_25km.tif')

# Directories
root <- '//CATALOGUE/WFP_ClimateRiskPr1'
aux  <- paste0(root,'/agroclimExtremes/agex_raw_data/agex_phenology')
out  <- paste0(root,'/agroclimExtremes/agex_raw_data'); dir.create(out, F, T)

# Mapping dekads to day-of-year units
yrs <- 1979:2023
doy <- lapply(X = yrs, FUN = function(yr){
  ## Select the median date within 10-days window
  # dys <- seq(from = as.Date(paste0(yr,'-01-01')), to = as.Date(paste0(yr,'-12-31')), by = 'day')
  # dks <- split(dys, ceiling(seq_along(dys)/10))
  # dks <- dks[1:36]
  # dks <- purrr::map(.x = dks, .f = function(x){as.character(median(x))}) |> unlist() |> as.Date()
  ## Select the initial date within 10-days window
  dks <- seq(from = as.Date(paste0(yr,'-01-01')), to = as.Date(paste0(yr,'-12-31')), length.out = 36)
  doy <- data.frame(dkd = lubridate::yday(dks))
  return(doy)
}) |> dplyr::bind_cols()
doy <- apply(X = doy, MARGIN = 1, FUN = function(x){round(mean(x))})
dks_map <- data.frame(day = doy, dekad = 1:36); rm(yrs, doy)
dks_map$day[dks_map$dekad == 36] <- 1

# Load number of growing seasons per site
n_seasons <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_phenology/phenonseasons_v03.tif')); gc(T)
n_seasons <- terra::resample(x = n_seasons, y = tmp_25km, method = 'near', threads = T); gc(T) # Growing seasons per site at 25 km resolution

# ------------------------------------------ #
# One season - Season 1 pre-processing
# ------------------------------------------ #
# Start date
if(!file.exists(paste0(aux,'/agex_one_s1_ini_raw_10km.tif'))){
  s1_ini <- terra::rast(paste0(aux,'/phenos1_v03.tif'))
  s1_ini[s1_ini > 250] <- NA # Remove missing data
  s1_ini_10km <- terra::resample(x = s1_ini, y = tmp_10km, method = 'near', threads = T)
  s1_ini_10km[n_seasons == 2] <- NA
  terra::writeRaster(x = s1_ini_10km, filename = paste0(aux,'/agex_one_s1_ini_raw_10km.tif'), overwrite = T)
  rm(s1_ini); gc(T)
}
s1_ini_10km <- terra::rast(paste0(aux,'/agex_one_s1_ini_raw_10km.tif'))
s1_ini_10km[s1_ini_10km > 36 & s1_ini_10km <= 72] <- s1_ini_10km[s1_ini_10km > 36 & s1_ini_10km <= 72] - 36 # Circularity
if(any(terra::values(s1_ini_10km) > 72, na.rm = T)){
  s1_ini_10km[s1_ini_10km > 72] <- s1_ini_10km[s1_ini_10km > 72] - 72
} # Circularity
s1_ini_10km <- terra::subst(s1_ini_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_ini_10km, filename = paste0(out,'/agex_one_s1_ini_10km.tif'), overwrite = T)
rm(s1_ini_10km); gc(T)

# End date
if(!file.exists(paste0(aux,'/agex_one_s1_end_raw_10km.tif'))){
  s1_end <- terra::rast(paste0(aux,'/phenoe1_v03.tif'))
  s1_end[s1_end > 250] <- NA
  s1_end_10km <- terra::resample(x = s1_end, y = tmp_10km, method = 'near', threads = T)
  s1_end_10km[n_seasons == 2] <- NA
  terra::writeRaster(x = s1_end_10km, filename = paste0(aux,'/agex_one_s1_end_raw_10km.tif'), overwrite = T)
  rm(s1_end); gc(T)
}
s1_end_10km <- terra::rast(paste0(aux,'/agex_one_s1_end_raw_10km.tif'))
s1_end_10km[s1_end_10km > 36 & s1_end_10km <= 72] <- s1_end_10km[s1_end_10km > 36 & s1_end_10km <= 72] - 36
if(any(terra::values(s1_end_10km) > 72, na.rm = T)){
  s1_end_10km[s1_end_10km > 72] <- s1_end_10km[s1_end_10km > 72] - 72
}
s1_end_10km <- terra::subst(s1_end_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_end_10km, filename = paste0(out,'/agex_one_s1_end_10km.tif'), overwrite = T)
rm(s1_end_10km); gc(T)

# ------------------------------------------ #
# Two seasons - Season 1 pre-processing
# ------------------------------------------ #
# Start date
if(!file.exists(paste0(aux,'/agex_two_s1_ini_raw_10km.tif'))){
  s1_ini <- terra::rast(paste0(aux,'/phenos1_v03.tif'))
  s1_ini[s1_ini > 250] <- NA # Remove missing data
  s1_ini_10km <- terra::resample(x = s1_ini, y = tmp_10km, method = 'near', threads = T)
  s1_ini_10km[n_seasons == 1] <- NA
  terra::writeRaster(x = s1_ini_10km, filename = paste0(aux,'/agex_two_s1_ini_raw_10km.tif'), overwrite = T)
  rm(s1_ini); gc(T)
}
s1_ini_10km <- terra::rast(paste0(aux,'/agex_two_s1_ini_raw_10km.tif'))
s1_ini_10km[s1_ini_10km > 36 & s1_ini_10km <= 72] <- s1_ini_10km[s1_ini_10km > 36 & s1_ini_10km <= 72] - 36 # Circularity
if(any(terra::values(s1_ini_10km) > 72, na.rm = T)){
  s1_ini_10km[s1_ini_10km > 72] <- s1_ini_10km[s1_ini_10km > 72] - 72
} # Circularity
s1_ini_10km <- terra::subst(s1_ini_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_ini_10km, filename = paste0(out,'/agex_two_s1_ini_10km.tif'), overwrite = T)
rm(s1_ini_10km); gc(T)

# End date
if(!file.exists(paste0(aux,'/agex_two_s1_end_raw_10km.tif'))){
  s1_end <- terra::rast(paste0(aux,'/phenoe1_v03.tif'))
  s1_end[s1_end > 250] <- NA
  s1_end_10km <- terra::resample(x = s1_end, y = tmp_10km, method = 'near', threads = T)
  s1_end_10km[n_seasons == 1] <- NA
  terra::writeRaster(x = s1_end_10km, filename = paste0(aux,'/agex_two_s1_end_raw_10km.tif'), overwrite = T)
  rm(s1_end); gc(T)
}
s1_end_10km <- terra::rast(paste0(aux,'/agex_two_s1_end_raw_10km.tif'))
s1_end_10km[s1_end_10km > 36 & s1_end_10km <= 72] <- s1_end_10km[s1_end_10km > 36 & s1_end_10km <= 72] - 36
if(any(terra::values(s1_end_10km) > 72, na.rm = T)){
  s1_end_10km[s1_end_10km > 72] <- s1_end_10km[s1_end_10km > 72] - 72
}
s1_end_10km <- terra::subst(s1_end_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_end_10km, filename = paste0(out,'/agex_two_s1_end_10km.tif'), overwrite = T)
rm(s1_end_10km); gc(T)

# ------------------------------------------ #
# Two seasons - Season 2 pre-processing
# ------------------------------------------ #
# Start date
if(!file.exists(paste0(aux,'/agex_two_s2_ini_raw_10km.tif'))){
  s2_ini <- terra::rast(paste0(aux,'/phenos2_v03.tif'))
  s2_ini[s2_ini > 250] <- NA # Remove missing data
  s2_ini_10km <- terra::resample(x = s2_ini, y = tmp_10km, method = 'near', threads = T)
  s2_ini_10km[n_seasons == 1] <- NA
  terra::writeRaster(x = s2_ini_10km, filename = paste0(aux,'/agex_two_s2_ini_raw_10km.tif'), overwrite = T)
  rm(s2_ini); gc(T)
}
s2_ini_10km <- terra::rast(paste0(aux,'/agex_two_s2_ini_raw_10km.tif'))
s2_ini_10km[s2_ini_10km > 36 & s2_ini_10km <= 72] <- s2_ini_10km[s2_ini_10km > 36 & s2_ini_10km <= 72] - 36 # Circularity
if(any(terra::values(s2_ini_10km) > 72, na.rm = T)){
  s2_ini_10km[s2_ini_10km > 72] <- s2_ini_10km[s2_ini_10km > 72] - 72
} # Circularity
s2_ini_10km <- terra::subst(s2_ini_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s2_ini_10km, filename = paste0(out,'/agex_two_s2_ini_10km.tif'), overwrite = T)
rm(s2_ini_10km); gc(T)

# End date
if(!file.exists(paste0(aux,'/agex_two_s2_end_raw_10km.tif'))){
  s2_end <- terra::rast(paste0(aux,'/phenoe2_v03.tif'))
  s2_end[s2_end > 250] <- NA
  s2_end_10km <- terra::resample(x = s2_end, y = tmp_10km, method = 'near', threads = T)
  s2_end_10km[n_seasons == 1] <- NA
  terra::writeRaster(x = s2_end_10km, filename = paste0(aux,'/agex_two_s2_end_raw_10km.tif'), overwrite = T)
  rm(s2_end); gc(T)
}
s2_end_10km <- terra::rast(paste0(aux,'/agex_two_s2_end_raw_10km.tif'))
s2_end_10km[s2_end_10km > 36 & s2_end_10km <= 72] <- s2_end_10km[s2_end_10km > 36 & s2_end_10km <= 72] - 36
if(any(terra::values(s2_end_10km) > 72, na.rm = T)){
  s2_end_10km[s2_end_10km > 72] <- s2_end_10km[s2_end_10km > 72] - 72
}
s2_end_10km <- terra::subst(s2_end_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s2_end_10km, filename = paste0(out,'/agex_two_s2_end_10km.tif'), overwrite = T)
rm(s2_end_10km); gc(T)

## ------------------------------------------ ##
## Filtering of growing seasons within croplands
## ------------------------------------------ ##

# One season
s1_ini_10km <- terra::rast(paste0(out,'/agex_one_s1_ini_10km.tif'))
croplnds    <- terra::rast(paste0(out,'/agex_croplands_mapspam_10km.tif'))
s1_ini_10km_msk <- terra::mask(x = s1_ini_10km, mask = croplnds)
terra::writeRaster(x = s1_ini_10km_msk, filename = paste0(out,'/agex_one_s1_ini_10km_croplands.tif'))
s1_end_10km <- terra::rast(paste0(out,'/agex_one_s1_end_10km.tif'))
s1_end_10km_msk <- terra::mask(x = s1_end_10km, mask = croplnds)
terra::writeRaster(x = s1_end_10km_msk, filename = paste0(out,'/agex_one_s1_end_10km_croplands.tif'))
rm(s1_ini_10km, s1_ini_10km_msk, s1_end_10km, s1_end_10km_msk); gc(T)

# Two seasons - Season 1
s1_ini_10km <- terra::rast(paste0(out,'/agex_two_s1_ini_10km.tif'))
s1_ini_10km_msk <- terra::mask(x = s1_ini_10km, mask = croplnds)
terra::writeRaster(x = s1_ini_10km_msk, filename = paste0(out,'/agex_two_s1_ini_10km_croplands.tif'))
s1_end_10km <- terra::rast(paste0(out,'/agex_two_s1_end_10km.tif'))
s1_end_10km_msk <- terra::mask(x = s1_end_10km, mask = croplnds)
terra::writeRaster(x = s1_end_10km_msk, filename = paste0(out,'/agex_two_s1_end_10km_croplands.tif'))
rm(s1_ini_10km, s1_ini_10km_msk, s1_end_10km, s1_end_10km_msk); gc(T)
# Two seasons - Season 2
s2_ini_10km <- terra::rast(paste0(out,'/agex_two_s2_ini_10km.tif'))
s2_ini_10km_msk <- terra::mask(x = s2_ini_10km, mask = croplnds)
terra::writeRaster(x = s2_ini_10km_msk, filename = paste0(out,'/agex_two_s2_ini_10km_croplands.tif'))
s2_end_10km <- terra::rast(paste0(out,'/agex_two_s2_end_10km.tif'))
s2_end_10km_msk <- terra::mask(x = s2_end_10km, mask = croplnds)
terra::writeRaster(x = s2_end_10km_msk, filename = paste0(out,'/agex_two_s2_end_10km_croplands.tif'))
rm(s2_ini_10km, s2_ini_10km_msk, s2_end_10km, s2_end_10km_msk); gc(T)

## ------------------------------------------ ##
## Growing seasons at 25 km
## ------------------------------------------ ##

## Number of growing seasons at 25 km
ngs <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_phenology/phenonseasons_v03.tif'))
tmp <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5_25km.tif')
ngs <- terra::resample(x = ngs, y = tmp, method = 'near')
terra::writeRaster(x = ngs, filename = paste0(root,'/agroclimExtremes/agex_raw_data/agex_nseasons_25km.tif'), overwrite = T)
