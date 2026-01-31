# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Resampling monthly SPEI-6 from 0.1° to 0.25°
# By: Harold Achicanoy
# WUR & ABC
# Created in March 2024
# Modified in January 2026
# --------------------------------------------------------------- #

# R options and user-defined functions
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,lubridate))

# Root directory
root <- '//CATALOGUE/AgroclimExtremes'

# Raster template at 25 km
tmp <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5_25km.tif')

# Resampling function
resampling_monthly_indices <- function(index = 'spei-6'){
  
  cat(paste0('>>> Resampling index: ',toupper(index),'... \n'))
  # List files
  fls <- list.files(path = paste0(root,'/agex_raw_data/monthly_',index,'_10km'), pattern = paste0('^',index,'*.*.tif$'), full.names = T)
  # Load indices at 10 km
  idx <- terra::rast(fls)
  names(idx) <- gsub('.tif','',basename(fls))
  # Croplands mask
  crp <- terra::rast(paste0(root,'/agex_raw_data/agex_croplands_foods_mapspam_25km.tif'))
  # Resample indices at 25 km
  idx_25km <- terra::resample(x = idx, y = tmp, method = 'max', threads = T)
  idx_25km <- terra::mask(x = idx_25km, mask = crp)
  outfiles <- paste0(root,'/agex_raw_data/monthly_',index,'_25km/',basename(fls))
  dir.create(dirname(outfiles[1]),F,T)
  terra::writeRaster(x = idx_25km, filename = outfiles, overwrite = T)
  
}
resampling_monthly_indices(index = 'spei-6')
