# ------------------------------------------ #
# Monthly SPEI from 0.1 to 0.25 degrees res
# By: Harold Achicanoy
# WUR & ABC
# Mar. 2024
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,lubridate))

root <- '//CATALOGUE/WFP_ClimateRiskPr1'

inp_path <- paste0(root,'/agroclimExtremes/agex_raw_data/monthly_spei_10km')
out_path <- paste0(root,'/agroclimExtremes/agex_raw_data/monthly_spei_25km')

dir.create(out_path, F, T)

# Template rasters
tmp_25km <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5_25km.tif')
cropland <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_croplands_foods_mapspam_25km.tif'))

fls <- list.files(path = inp_path, pattern = '.tif$', full.names = T)
nms <- basename(fls)

lapply(1:length(fls), function(i){
  r <- terra::rast(fls[i])
  r <- terra::resample(x = r, y = tmp_25km, method = 'max', threads = T)
  r <- terra::mask(x = r, mask = cropland)
  terra::writeRaster(x = r, filename = paste0(out_path,'/',nms[i]), overwrite = T)
})
