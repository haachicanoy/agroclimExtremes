# ------------------------------------------ #
# Prepare study areas
# By: Harold Achicanoy
# WUR & ABC
# Dec. 2023
# ------------------------------------------ #

options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,tidyverse,geodata))

# Define directories
inp_dir <- 'D:/Data/Maps'
out_dir <- 'D:/OneDrive - CGIAR/PhD/papers/paper1'

# Cropland areas from CROPGRIDS
crops_dir <- paste0(inp_dir,'/cropgrids/NC_maps') # Directory
crops_fls <- list.files(path = crops_dir, pattern = '.nc', full.names = T)
crops_fls <- crops_fls[-1] # File names

crops_lnd <- purrr::map(.x = crops_fls, .f = function(x){ # Data
  r <- terra::rast(x)[[1]]; return(r)
}) |> terra::rast() |> sum()

crops_lnd[crops_lnd <= 0] <- NA
crops_lnd[!is.na(crops_lnd)] <- 1

terra::writeRaster(x = crops_lnd, filename = paste0(out_dir,'/data/cropgrids_harvested_area.tif'))

# Cropland areas from MAPSPAM 2010
crops_dir <- paste0(inp_dir,'/spam2010') # Directory
crops_fls <- list.files(path = crops_dir, pattern = '_A.tif', full.names = T)

crops_lnd <- terra::rast(crops_fls) |> sum(na.rm = T)
crops_lnd[crops_lnd <= 0] <- NA
crops_lnd[!is.na(crops_lnd)] <- 1

terra::writeRaster(x = crops_lnd, filename = paste0(out_dir,'/data/mapspam_harvested_area.tif'), overwrite = T)

# Mask of water bodies and protected areas
wbd <- terra::rast(paste0(inp_dir,'/water_bodies/glwd_3')) # Water bodies
pra <- terra::rast(list.files(path = paste0(inp_dir,'/protected_areas'), pattern = '.tif$', full.names = T)) # Protected areas
pra <- sum(pra, na.rm = T)

msk <- c(wbd, pra); rm(wbd, pra)
msk <- sum(msk, na.rm = T)

msk[msk > 0] <- 1

terra::writeRaster(x = msk,
                   filename = paste0(out_dir,'/data/mask_water_protectedareas.tif'),
                   gdal = c('COMPRESS=DEFLATE', 'TFW=YES'))

# Masking MAPSPAM by water bodies and protected areas
mapspam <- terra::rast(paste0(out_dir,'/data/mapspam_harvested_area.tif'))
msk_10 <- terra::resample(x = msk, y = mapspam, method = 'near')
mapspam_mskd <- terra::mask(x = mapspam, mask = msk_10, inverse = T)
terra::writeRaster(x = mapspam_mskd,
                   filename = paste0(out_dir,'/data/mapspam_harvested_area_masked.tif'))

geodata::cropland(source = 'WorldCover', path = paste0(out_dir,'/data'))
geodata::cropland(source = 'GLAD', year = 2019, path = paste0(out_dir,'/data'))
geodata::cropland(source = 'QED', path = paste0(out_dir,'/data'))

# Cattle areas
lvstc_lnd <- terra::rast('D:/OneDrive - CGIAR/African_Crisis_Observatory/data/_global/livestock/cattle/5_Ct_2010_Da.tif')
lvstc_lnd <- lvstc_lnd >= 100
plot(lvstc_lnd)
plot(mapspam)
