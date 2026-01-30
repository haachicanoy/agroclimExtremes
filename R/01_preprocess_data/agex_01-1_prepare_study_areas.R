# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Prepare agricultural production areas
# By: Harold Achicanoy
# WUR & ABC
# Created in December 2023
# Modified in January 2026
# --------------------------------------------------------------- #

# R options and user-defined functions
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra))
list.files2 <- Vectorize(FUN = list.files, vectorize.args = 'path')
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')

# Template rasters at 0.1° and 0.25°
tmp_10km <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5.tif')
tmp_25km <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5_25km.tif')

# Define directories
root <- '//CATALOGUE/AgroclimExtremes'
inp_dir <- paste0(root,'/agex_raw_data/croplands')

## Cropland areas from MapSPAM 2010
crops_dir <- paste0(inp_dir,'/spam2010') # Directory
crops_fls <- list.files(path = crops_dir, pattern = '_A.tif', full.names = T)

# Crops classification
grp <- read.csv('https://raw.githubusercontent.com/wri/MAPSPAM/master/metadata_tables/4-Methodology-Crops-of-SPAM-2005-2015-02-26.csv')
grp <- grp[,c('SPAM.short.name','GROUP')]
foods <- grp[!(grp$GROUP %in% c('fibres','stimulant')),'SPAM.short.name']
crops_fls <- crops_fls[grep2(pattern = toupper(foods), x = crops_fls)]

crops_lnd <- terra::rast(crops_fls) |> sum(na.rm = T)
crops_lnd[crops_lnd <= 0] <- NA
crops_lnd[!is.na(crops_lnd)] <- 1

# Resampling MapSPAM into template rasters resolution
crops_lnd_10km <- terra::resample(x = crops_lnd, y = tmp_10km, method = 'near')
terra::writeRaster(x = crops_lnd_10km, filename = paste0(root,'/agex_raw_data/agex_croplands_foods_mapspam_10km.tif'), overwrite = T)
crops_lnd_25km <- terra::resample(x = crops_lnd, y = tmp_25km, method = 'near')
terra::writeRaster(x = crops_lnd_25km, filename = paste0(root,'/agex_raw_data/agex_croplands_foods_mapspam_25km.tif'), overwrite = T)
