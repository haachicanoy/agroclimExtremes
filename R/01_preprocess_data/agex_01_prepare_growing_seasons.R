# ------------------------------------------ #
# Prepare growing seasons - croplands
# By: Harold Achicanoy
# WUR & ABC
# Dec. 2023
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,fields,spatstat,future,furrr))

# Template rasters
# tmp_10km <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5.tif')
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
if(!file.exists(paste0(out,'/agex_nseasons_25km_corrected.tif'))){
  n_seasons <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_phenology/phenonseasons_v03.tif')); gc(T)
  n_seasons <- terra::resample(x = n_seasons, y = tmp_25km, method = 'near', threads = T); gc(T) # Growing seasons per site at 25 km resolution
  croplands <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_croplands_foods_mapspam_25km.tif'))
  n_seasons <- terra::mask(x = n_seasons, mask = croplands)
  terra::writeRaster(x = n_seasons, filename = paste0(out,'/agex_nseasons_25km.tif'), overwrite = T)
  
  # Split seasons
  n_seasons_1 <- n_seasons # Season 1
  n_seasons_1[n_seasons_1 != 1] <- NA
  n_seasons_1_dfm <- terra::as.data.frame(x = n_seasons_1, xy = T, cell = T, na.rm = T) # Coordinates
  
  n_seasons_2 <- n_seasons # Season 2
  n_seasons_2[n_seasons_2 != 2] <- NA
  n_seasons_2_dfm <- terra::as.data.frame(x = n_seasons_2, xy = T, cell = T, na.rm = T)
  
  r_owin <- as.owin(unlist(as.list(terra::ext(n_seasons))))
  
  # Get probabilities from Season 1 coordinates
  p1 <- ppp(x = n_seasons_1_dfm$x, y = n_seasons_1_dfm$y, window = r_owin)
  ds1 <- density(p1)
  ds_r1 <- terra::rast(ds1)
  ds_r1 <- terra::resample(x = ds_r1, y = n_seasons_1)
  ds_r1 <- terra::mask(x = ds_r1, mask = n_seasons_1)
  ds_r1 <- ds_r1/sum(terra::values(ds_r1, na.rm = T))
  
  # Get probabilities from Season 2 coordinates
  p2 <- ppp(x = n_seasons_2_dfm$x, y = n_seasons_2_dfm$y, window = r_owin)
  ds2 <- density(p2)
  ds_r2 <- terra::rast(ds2)
  ds_r2 <- terra::resample(x = ds_r2, y = n_seasons_2)
  ds_r2 <- terra::mask(x = ds_r2, mask = n_seasons_2)
  ds_r2 <- ds_r2/sum(terra::values(ds_r2, na.rm = T))
  
  names(n_seasons_1_dfm)[ncol(n_seasons_1_dfm)] <- 'nseasons'
  names(n_seasons_2_dfm)[ncol(n_seasons_2_dfm)] <- 'nseasons'
  
  n_seasons_dfm <- rbind(n_seasons_1_dfm, n_seasons_2_dfm)
  
  probs <- rbind(terra::as.data.frame(x = ds_r1, xy = T, cell = T, na.rm = T),
                 terra::as.data.frame(x = ds_r2, xy = T, cell = T, na.rm = T))
  names(probs)[ncol(probs)] <- 'probability'
  probs$probability <- probs$probability/2
  
  n_seasons_dfm <- dplyr::left_join(x = n_seasons_dfm, y = probs[,c('cell','probability')], by = 'cell')
  rm(probs, ds_r1, ds_r2, ds1, ds2, p1, p2, r_owin, n_seasons_1_dfm, n_seasons_2_dfm, n_seasons_1, n_seasons_2)
  
  # Visualize outliers
  hist(n_seasons_dfm$probability[n_seasons_dfm$nseasons == 2])
  thr <- 0.000017
  abline(v = thr, col = 2, lty = 2)
  plot(n_seasons)
  points(n_seasons_dfm[n_seasons_dfm$nseasons == 2 & n_seasons_dfm$probability < thr,c('x','y')],
         col = 'red', pch = 20)
  
  n_seasons_dfm2 <- n_seasons_dfm
  n_seasons_dfm2$nseasons[n_seasons_dfm2$nseasons == 2 & n_seasons_dfm2$probability < thr] <- 1
  
  aux <- n_seasons
  terra::values(aux) <- NA
  aux[n_seasons_dfm2$cell] <- n_seasons_dfm2$nseasons
  
  terra::writeRaster(x = aux, filename = paste0(out,'/agex_nseasons_25km_corrected.tif'), overwrite = T)
  saveRDS(object = n_seasons_dfm, file = paste0(out,'/agex_nseasons_25km_dfm.rds'))
  rm(n_seasons_dfm2, n_seasons_dfm, aux, croplands, thr)
  
  # hist(n_seasons_dfm$probability[n_seasons_dfm$nseasons == 1])
  # thr <- 0.0000011
  # abline(v = thr, col = 2, lty = 2)
  # plot(-1*(n_seasons-2))
  # points(n_seasons_dfm[n_seasons_dfm$nseasons == 1 & n_seasons_dfm$probability < thr,c('x','y')],
  #        col = 'red', pch = 20)
} else {
  n_seasons <- terra::rast(paste0(out,'/agex_nseasons_25km_corrected.tif'))
  n_seasons_dfm <- readRDS(file = paste0(out,'/agex_nseasons_25km_dfm.rds'))
  thr <- 0.000017
}

# clls <- n_seasons_dfm$cell[n_seasons_dfm$nseasons == 2 & n_seasons_dfm$probability < thr]
# 
# adj_clls <- terra::adjacent(x = n_seasons, cells = clls, directions = 'queen', pairs = F, symmetrical = T)
# terra::extract(x = n_seasons, y = terra::xyFromCell(n_seasons, cell = adj_clls[6,]))

# ------------------------------------------ #
# One season - Season 1 pre-processing
# ------------------------------------------ #
# Start date
if(!file.exists(paste0(aux,'/agex_one_s1_ini_raw_25km.tif'))){
  s1_ini <- terra::rast(paste0(aux,'/phenos1_v03.tif'))
  s1_ini[s1_ini > 250] <- NA # Remove missing data
  s1_ini_25km <- terra::resample(x = s1_ini, y = tmp_25km, method = 'near', threads = T)
  s1_ini_25km[n_seasons == 2] <- NA
  terra::writeRaster(x = s1_ini_25km, filename = paste0(aux,'/agex_one_s1_ini_raw_25km.tif'), overwrite = T)
  rm(s1_ini); gc(T)
}
s1_ini_25km <- terra::rast(paste0(aux,'/agex_one_s1_ini_raw_25km.tif'))
s1_ini_25km[s1_ini_25km > 36 & s1_ini_25km <= 72] <- s1_ini_25km[s1_ini_25km > 36 & s1_ini_25km <= 72] - 36 # Circularity
if(any(terra::values(s1_ini_25km) > 72, na.rm = T)){
  s1_ini_25km[s1_ini_25km > 72] <- s1_ini_25km[s1_ini_25km > 72] - 72
} # Circularity
s1_ini_25km <- terra::subst(s1_ini_25km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_ini_25km, filename = paste0(out,'/agex_one_s1_ini_25km.tif'), overwrite = T)
rm(s1_ini_25km); gc(T)

# End date
if(!file.exists(paste0(aux,'/agex_one_s1_end_raw_25km.tif'))){
  s1_end <- terra::rast(paste0(aux,'/phenoe1_v03.tif'))
  s1_end[s1_end > 250] <- NA
  s1_end_25km <- terra::resample(x = s1_end, y = tmp_25km, method = 'near', threads = T)
  s1_end_25km[n_seasons == 2] <- NA
  terra::writeRaster(x = s1_end_25km, filename = paste0(aux,'/agex_one_s1_end_raw_25km.tif'), overwrite = T)
  rm(s1_end); gc(T)
}
s1_end_25km <- terra::rast(paste0(aux,'/agex_one_s1_end_raw_25km.tif'))
s1_end_25km[s1_end_25km > 36 & s1_end_25km <= 72] <- s1_end_25km[s1_end_25km > 36 & s1_end_25km <= 72] - 36
if(any(terra::values(s1_end_25km) > 72, na.rm = T)){
  s1_end_25km[s1_end_25km > 72] <- s1_end_25km[s1_end_25km > 72] - 72
}
s1_end_25km <- terra::subst(s1_end_25km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_end_25km, filename = paste0(out,'/agex_one_s1_end_25km.tif'), overwrite = T)
rm(s1_end_25km); gc(T)

# ------------------------------------------ #
# Two seasons - Season 1 pre-processing
# ------------------------------------------ #
# Start date
if(!file.exists(paste0(aux,'/agex_two_s1_ini_raw_25km.tif'))){
  s1_ini <- terra::rast(paste0(aux,'/phenos1_v03.tif'))
  s1_ini[s1_ini > 250] <- NA # Remove missing data
  s1_ini_25km <- terra::resample(x = s1_ini, y = tmp_25km, method = 'near', threads = T)
  s1_ini_25km[n_seasons == 1] <- NA
  terra::writeRaster(x = s1_ini_25km, filename = paste0(aux,'/agex_two_s1_ini_raw_25km.tif'), overwrite = T)
  rm(s1_ini); gc(T)
}
s1_ini_25km <- terra::rast(paste0(aux,'/agex_two_s1_ini_raw_25km.tif'))
s1_ini_25km[s1_ini_25km > 36 & s1_ini_25km <= 72] <- s1_ini_25km[s1_ini_25km > 36 & s1_ini_25km <= 72] - 36 # Circularity
if(any(terra::values(s1_ini_25km) > 72, na.rm = T)){
  s1_ini_25km[s1_ini_25km > 72] <- s1_ini_25km[s1_ini_25km > 72] - 72
} # Circularity
s1_ini_25km <- terra::subst(s1_ini_25km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_ini_25km, filename = paste0(out,'/agex_two_s1_ini_25km.tif'), overwrite = T)
rm(s1_ini_25km); gc(T)

# End date
if(!file.exists(paste0(aux,'/agex_two_s1_end_raw_25km.tif'))){
  s1_end <- terra::rast(paste0(aux,'/phenoe1_v03.tif'))
  s1_end[s1_end > 250] <- NA
  s1_end_25km <- terra::resample(x = s1_end, y = tmp_25km, method = 'near', threads = T)
  s1_end_25km[n_seasons == 1] <- NA
  terra::writeRaster(x = s1_end_25km, filename = paste0(aux,'/agex_two_s1_end_raw_25km.tif'), overwrite = T)
  rm(s1_end); gc(T)
}
s1_end_25km <- terra::rast(paste0(aux,'/agex_two_s1_end_raw_25km.tif'))
s1_end_25km[s1_end_25km > 36 & s1_end_25km <= 72] <- s1_end_25km[s1_end_25km > 36 & s1_end_25km <= 72] - 36
if(any(terra::values(s1_end_25km) > 72, na.rm = T)){
  s1_end_25km[s1_end_25km > 72] <- s1_end_25km[s1_end_25km > 72] - 72
}
s1_end_25km <- terra::subst(s1_end_25km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_end_25km, filename = paste0(out,'/agex_two_s1_end_25km.tif'), overwrite = T)
rm(s1_end_25km); gc(T)

# ------------------------------------------ #
# Two seasons - Season 2 pre-processing
# ------------------------------------------ #
# Start date
if(!file.exists(paste0(aux,'/agex_two_s2_ini_raw_25km.tif'))){
  s2_ini <- terra::rast(paste0(aux,'/phenos2_v03.tif'))
  s2_ini[s2_ini > 250] <- NA # Remove missing data
  s2_ini_25km <- terra::resample(x = s2_ini, y = tmp_25km, method = 'near', threads = T)
  s2_ini_25km[n_seasons == 1] <- NA
  terra::writeRaster(x = s2_ini_25km, filename = paste0(aux,'/agex_two_s2_ini_raw_25km.tif'), overwrite = T)
  rm(s2_ini); gc(T)
}
s2_ini_25km <- terra::rast(paste0(aux,'/agex_two_s2_ini_raw_25km.tif'))
s2_ini_25km[s2_ini_25km > 36 & s2_ini_25km <= 72] <- s2_ini_25km[s2_ini_25km > 36 & s2_ini_25km <= 72] - 36 # Circularity
if(any(terra::values(s2_ini_25km) > 72, na.rm = T)){
  s2_ini_25km[s2_ini_25km > 72] <- s2_ini_25km[s2_ini_25km > 72] - 72
} # Circularity
s2_ini_25km <- terra::subst(s2_ini_25km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s2_ini_25km, filename = paste0(out,'/agex_two_s2_ini_25km.tif'), overwrite = T)
rm(s2_ini_25km); gc(T)

# End date
if(!file.exists(paste0(aux,'/agex_two_s2_end_raw_25km.tif'))){
  s2_end <- terra::rast(paste0(aux,'/phenoe2_v03.tif'))
  s2_end[s2_end > 250] <- NA
  s2_end_25km <- terra::resample(x = s2_end, y = tmp_25km, method = 'near', threads = T)
  s2_end_25km[n_seasons == 1] <- NA
  terra::writeRaster(x = s2_end_25km, filename = paste0(aux,'/agex_two_s2_end_raw_25km.tif'), overwrite = T)
  rm(s2_end); gc(T)
}
s2_end_25km <- terra::rast(paste0(aux,'/agex_two_s2_end_raw_25km.tif'))
s2_end_25km[s2_end_25km > 36 & s2_end_25km <= 72] <- s2_end_25km[s2_end_25km > 36 & s2_end_25km <= 72] - 36
if(any(terra::values(s2_end_25km) > 72, na.rm = T)){
  s2_end_25km[s2_end_25km > 72] <- s2_end_25km[s2_end_25km > 72] - 72
}
s2_end_25km <- terra::subst(s2_end_25km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s2_end_25km, filename = paste0(out,'/agex_two_s2_end_25km.tif'), overwrite = T)
rm(s2_end_25km); gc(T)

## ------------------------------------------ ##
## Filtering of growing seasons within croplands
## ------------------------------------------ ##

# One season
s1_ini_25km <- terra::rast(paste0(out,'/agex_one_s1_ini_25km.tif'))
croplnds    <- terra::rast(paste0(out,'/agex_croplands_foods_mapspam_25km.tif'))
s1_ini_25km_msk <- terra::mask(x = s1_ini_25km, mask = croplnds)
terra::writeRaster(x = s1_ini_25km_msk, filename = paste0(out,'/agex_one_s1_ini_25km_croplands.tif'))
s1_end_25km <- terra::rast(paste0(out,'/agex_one_s1_end_25km.tif'))
s1_end_25km_msk <- terra::mask(x = s1_end_25km, mask = croplnds)
terra::writeRaster(x = s1_end_25km_msk, filename = paste0(out,'/agex_one_s1_end_25km_croplands.tif'))
rm(s1_ini_25km, s1_ini_25km_msk, s1_end_25km, s1_end_25km_msk); gc(T)

# Two seasons - Season 1
s1_ini_25km <- terra::rast(paste0(out,'/agex_two_s1_ini_25km.tif'))
s1_ini_25km_msk <- terra::mask(x = s1_ini_25km, mask = croplnds)
terra::writeRaster(x = s1_ini_25km_msk, filename = paste0(out,'/agex_two_s1_ini_25km_croplands.tif'))
s1_end_25km <- terra::rast(paste0(out,'/agex_two_s1_end_25km.tif'))
s1_end_25km_msk <- terra::mask(x = s1_end_25km, mask = croplnds)
terra::writeRaster(x = s1_end_25km_msk, filename = paste0(out,'/agex_two_s1_end_25km_croplands.tif'))
rm(s1_ini_25km, s1_ini_25km_msk, s1_end_25km, s1_end_25km_msk); gc(T)
# Two seasons - Season 2
s2_ini_25km <- terra::rast(paste0(out,'/agex_two_s2_ini_25km.tif'))
s2_ini_25km_msk <- terra::mask(x = s2_ini_25km, mask = croplnds)
terra::writeRaster(x = s2_ini_25km_msk, filename = paste0(out,'/agex_two_s2_ini_25km_croplands.tif'))
s2_end_25km <- terra::rast(paste0(out,'/agex_two_s2_end_25km.tif'))
s2_end_25km_msk <- terra::mask(x = s2_end_25km, mask = croplnds)
terra::writeRaster(x = s2_end_25km_msk, filename = paste0(out,'/agex_two_s2_end_25km_croplands.tif'))
rm(s2_ini_25km, s2_ini_25km_msk, s2_end_25km, s2_end_25km_msk); gc(T)
