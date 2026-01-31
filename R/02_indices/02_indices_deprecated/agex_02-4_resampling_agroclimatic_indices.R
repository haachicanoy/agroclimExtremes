# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Resampling monthly SPEI-6 at 0.25°
# By: Harold Achicanoy
# WUR & ABC
# Created in February 2024
# Modified in January 2026
# --------------------------------------------------------------- #

# R options and user-defined functions
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,lubridate))

root <- '//CATALOGUE/AgroclimExtremes'

# Raster template at 25 km
tmp <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5_25km.tif')

# Resampling function
resampling_indices <- function(index = 'spei-6'){
  
  cat(paste0('>>> Resampling index: ',index,'... \n'))
  cat(paste0('> One growing season\n'))
  # List files
  fls <- list.files(path = paste0(root,'/agex_indices/agex_',index,'/agex_',index,'_10km'), pattern = paste0('^one_s1_',index,'_[0-9]*.*.tif$'), full.names = T)
  # Load indices at 10 km
  idx <- terra::rast(fls)
  # Resample indices at 25 km
  idx_25km <- terra::resample(x = idx, y = tmp, method = 'max')
  outfile <- paste0(root,'/agex_indices/agex_',index,'/agex_',index,'_25km/one_s1_',index,'_25km.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = idx_25km, filename = outfile, overwrite = T)
  
  cat(paste0('> Two growing seasons: season 1.\n'))
  # List files
  fls <- list.files(path = paste0(root,'/agex_indices/agex_',index,'/agex_',index,'_10km'), pattern = paste0('^two_s1_',index,'_[0-9]*.*.tif$'), full.names = T)
  # Load indices at 10 km
  idx <- terra::rast(fls)
  # Resample indices at 25 km
  idx_25km <- terra::resample(x = idx, y = tmp, method = 'max')
  outfile <- paste0(root,'/agex_indices/agex_',index,'/agex_',index,'_25km/two_s1_',index,'_25km.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = idx_25km, filename = outfile, overwrite = T)
  
  cat(paste0('> Two growing seasons: season 2.\n'))
  # List files
  fls <- list.files(path = paste0(root,'/agex_indices/agex_',index,'/agex_',index,'_10km'), pattern = paste0('^two_s2_',index,'_[0-9]*.*.tif$'), full.names = T)
  # Load indices at 10 km
  idx <- terra::rast(fls)
  # Resample indices at 25 km
  idx_25km <- terra::resample(x = idx, y = tmp, method = 'max')
  outfile <- paste0(root,'/agex_indices/agex_',index,'/agex_',index,'_25km/two_s2_',index,'_25km.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = idx_25km, filename = outfile, overwrite = T)
  cat('Done!\n\n')
  
}

resampling_indices(index = 'spei-6') # SPEI-6
