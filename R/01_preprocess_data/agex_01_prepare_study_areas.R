# ------------------------------------------ #
# Prepare study areas
# By: Harold Achicanoy
# WUR & ABC
# Dec. 2023
# ------------------------------------------ #

options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,tidyverse,geodata,factoextra))

list.files2 <- Vectorize(FUN = list.files, vectorize.args = 'path')
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')

# Template rasters
tmp_10km <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5.tif')
tmp_25km <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5_25km.tif')

# Define directories
root <- '//CATALOGUE/WFP_ClimateRiskPr1'           # Server
inp_dir <- 'D:/Data/Maps'                          # Local
out_dir <- 'D:/OneDrive - CGIAR/PhD/papers/paper1' # local

# # Cropland areas from CROPGRIDS
# crops_dir <- paste0(inp_dir,'/cropgrids/NC_maps') # Directory
# crops_fls <- list.files(path = crops_dir, pattern = '.nc', full.names = T)
# crops_fls <- crops_fls[-1] # File names
# 
# crops_lnd <- purrr::map(.x = crops_fls, .f = function(x){ # Data
#   r <- terra::rast(x)[[1]]; return(r)
# }) |> terra::rast() |> sum()
# 
# crops_lnd[crops_lnd <= 0] <- NA
# crops_lnd[!is.na(crops_lnd)] <- 1
# 
# terra::writeRaster(x = crops_lnd, filename = paste0(out_dir,'/data/croplands_cropgrids.tif'))

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

terra::writeRaster(x = crops_lnd, filename = paste0(out_dir,'/data/agex_croplands_foods_mapspam.tif'), overwrite = T)

# Resampling MapSPAM into AgERA5 template resolution
crops_lnd_10km <- terra::resample(x = crops_lnd, y = tmp_10km, method = 'near')
terra::writeRaster(x = crops_lnd_10km, filename = paste0(root,'/agroclimExtremes/agex_raw_data/agex_croplands_foods_mapspam_10km.tif'), overwrite = T)
crops_lnd_25km <- terra::resample(x = crops_lnd, y = tmp_25km, method = 'near')
terra::writeRaster(x = crops_lnd_25km, filename = paste0(root,'/agroclimExtremes/agex_raw_data/agex_croplands_foods_mapspam_25km.tif'), overwrite = T)

# # Mask of water bodies and protected areas
# wbd <- terra::rast(paste0(inp_dir,'/water_bodies/glwd_3')) # Water bodies
# pra <- terra::rast(list.files(path = paste0(inp_dir,'/protected_areas'), pattern = '.tif$', full.names = T)) # Protected areas
# pra <- sum(pra, na.rm = T)
# 
# msk <- c(wbd, pra); rm(wbd, pra)
# msk <- sum(msk, na.rm = T)
# 
# msk[msk > 0] <- 1
# 
# terra::writeRaster(x = msk, filename = paste0(out_dir,'/data/mask_water_protectedareas.tif'))

# # Masking MAPSPAM by water bodies and protected areas
# mapspam <- terra::rast(paste0(out_dir,'/data/mapspam_harvested_area.tif'))
# msk_10 <- terra::resample(x = msk, y = mapspam, method = 'near')
# mapspam_mskd <- terra::mask(x = mapspam, mask = msk_10, inverse = T)
# terra::writeRaster(x = mapspam_mskd,
#                    filename = paste0(out_dir,'/data/mapspam_harvested_area_masked.tif'))
# 
# geodata::cropland(source = 'WorldCover', path = paste0(out_dir,'/data'))
# geodata::cropland(source = 'GLAD', year = 2019, path = paste0(out_dir,'/data'))
# geodata::cropland(source = 'QED', path = paste0(out_dir,'/data'))

# Global livestock production systems (GLPS). Not used
lvstc_sys <- terra::rast(paste0(inp_dir,'/livestock_systems/glps_gleam_61113_10km.tif'))
lvstc_sys[lvstc_sys == 15] <- NA # Removing unsuitable areas

# Animal areas
lvs_dir <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory/data/_global/livestock'
anmls   <- list.dirs(path = lvs_dir, full.names = F, recursive = F)
lvstc_fls <- list.files2(path = paste0(lvs_dir,'/',anmls), pattern = '_Da.tif$', full.names = T); rm(lvs_dir, anmls)
lvstc_cnt <- terra::rast(lvstc_fls)

# Computing Livestock Units
# LU: Livestock Unit
# LU = Buffaloes * 1 + Cattle * 1 + Chickens * ((0.007 + 0.014)/2)  + Ducks * 0.01 + Goats * 0.1 +
#      Horses * 0.8 + Pigs * ((0.5+0.027)/2) + Sheep * 0.1
# Source: https://ec.europa.eu/eurostat/statistics-explained/index.php?title=Glossary:Livestock_unit_(LSU)
lsu <- lvstc_cnt[[1]] + lvstc_cnt[[2]] + lvstc_cnt[[3]]*((0.007 + 0.014)/2) +
  lvstc_cnt[[4]]*0.01 + lvstc_cnt[[5]]*0.1 + lvstc_cnt[[6]]*0.8 +
  lvstc_cnt[[7]]*((0.5+0.027)/2) + lvstc_cnt[[8]]*0.1
rm(lvstc_fls, lvstc_cnt)

vls <- terra::values(lsu, na.rm = T) |> as.numeric()
hist(vls[vls > 0])
round(stats::quantile(x = vls[vls > 0], probs = seq(0,1,0.1)), 2)
plot(lsu > 0)
plot(lsu > 100)    # Selected threshold
plot(lsu > 145.96) # Median (distribution > 0)
rm(vls)

lsu <- (lsu > 100) * 1; lsu[lsu == 0] <- NA
# Resampling Livestock Units into AgERA5 template resolution
lvstc_lnd_10km <- terra::resample(x = lsu, y = tmp_10km, method = 'near')
terra::writeRaster(x = lvstc_lnd_10km, filename = paste0(root,'/agroclimExtremes/agex_raw_data/agex_livestockunits_FAO_10km.tif'), overwrite = T)

# # Replacing 0's with the minimum value of the distribution
# for(i in 1:(terra::nlyr(lvstc_cnt))){
#   vls <- as.numeric(terra::values(lvstc_cnt[[i]], na.rm = T))
#   vls_min <- min(vls[vls > 0])
#   lvstc_cnt[[i]][lvstc_cnt[[i]] == 0] <- vls_min; rm(vls_min)
# }
# lvstc_cnt <- log(lvstc_cnt) # Apply logarithm transformation to correct skewness
# lvstc_pca <- terra::prcomp(lvstc_cnt, center = T, scale = T) # PCA
# lvstc_pca
# 
# get_pca_var(lvstc_pca)$contrib
# 
# names(lvstc_cnt) <- paste0('X',names(lvstc_cnt))
# 
# lvstc_res <- terra::predict(lvstc_cnt, lvstc_pca)
# plot(lvstc_res[[1]] > 0)
# plot(lvstc_res[[2]] > 0)
# plot(lvstc_res[[6]] > 0)
