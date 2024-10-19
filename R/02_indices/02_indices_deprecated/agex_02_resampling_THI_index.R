## ------------------------------------------ ##
## Resampling THI index at 25 km
## By: Harold Achicanoy
## WUR & ABC
## Mar. 2024
## ------------------------------------------ ##

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,lubridate))

root <- '//CATALOGUE/AgroclimExtremes'

# Raster template at 25 km
tmp <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5_25km.tif')

# Resampling function
resampling_indices <- function(index = 'thi'){
  
  cat(paste0('>>> Resampling index: ',toupper(index),'... \n'))
  # List files
  fls <- list.files(path = paste0(root,'/agex_raw_data/monthly_',index,'_10km'), pattern = paste0('^',toupper(index),'_max-*.*.tif$'), full.names = T)
  # Load indices at 10 km
  idx <- terra::rast(fls)
  names(idx) <- gsub('.tif','',basename(fls))
  # Croplands mask
  crp <- terra::rast(paste0(root,'/agex_raw_data/agex_croplands_foods_mapspam_25km.tif'))
  # Resample indices at 25 km
  idx_25km <- terra::resample(x = idx, y = tmp, method = 'near', threads = T)
  idx_25km <- terra::mask(x = idx_25km, mask = crp)
  outfile <- paste0(root,'/agex_indices/agex_',index,'/agex_',index,'_25km/monthly_',index,'_25km.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = idx_25km, filename = outfile, overwrite = T)
  
  # Annual
  idx_annual_25km <- terra::tapp(x = idx_25km, index = gsub('-[0-9][0-9]','',gsub('THI_max-','',names(idx_25km))), fun = max)
  names(idx_annual_25km) <- paste0('THI_',1979:2023)
  outfile <- paste0(root,'/agex_indices/agex_',index,'/agex_',index,'_25km/annual_',index,'_25km.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = idx_annual_25km, filename = outfile, overwrite = T)
}
