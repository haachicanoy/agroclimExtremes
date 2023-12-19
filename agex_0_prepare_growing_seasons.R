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

# Directories
root <- '//CATALOGUE/WFP_ClimateRiskPr1/agroclimExtremes/data' # root <- 'D:/Data/Maps'
# root <- 'common_data/agroclimExtremes/data'
out  <- paste0(root,'/results'); dir.create(out, F, T)

# Mapping dekads to day-of-year units
yrs <- 1981:2023
doy <- lapply(X = yrs, FUN = function(yr){
  dys <- seq(from = as.Date(paste0(yr,'-01-01')), to = as.Date(paste0(yr,'-12-31')), by = 'day')
  dks <- split(dys, ceiling(seq_along(dys)/10))
  dks <- dks[1:36]
  dks <- purrr::map(.x = dks, .f = function(x){as.character(median(x))}) |> unlist() |> as.Date()
  doy <- data.frame(dkd = lubridate::yday(dks))
  return(doy)
}) |> dplyr::bind_cols()
doy <- apply(X = doy, MARGIN = 1, FUN = function(x){round(mean(x))})
dks_map <- data.frame(day = doy, dekad = 1:36); rm(yrs, doy)

# Load phenology data
n_seasons <- terra::rast(paste0(root,'/phenology/phenonseasons_v03.tif')); gc(TRUE)

# Template rasters
# tmp_05km <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_chirps.tif')
tmp_10km <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5.tif')

# ------------------------------------------ #
# One season - Season 1 pre-processing
# ------------------------------------------ #
# Start date
if(!file.exists(paste0(out,'/one_s1_ini_raw_10km.tif'))){
  s1_ini <- terra::rast(paste0(root,'/phenology/phenos1_v03.tif'))
  s1_ini[s1_ini > 250] <- NA # Remove missing data
  s1_ini[n_seasons == 2] <- NA
  s1_ini_10km <- terra::resample(x = s1_ini, y = tmp_10km, method = 'near', threads = T)
  rm(s1_ini); gc(TRUE)
  terra::writeRaster(x = s1_ini_10km, filename = paste0(out,'/one_s1_ini_raw_10km.tif'), overwrite = T)
}
s1_ini_10km <- terra::rast(paste0(out,'/one_s1_ini_raw_10km.tif'))
s1_ini_10km[s1_ini_10km > 36 & s1_ini_10km <= 72] <- s1_ini_10km[s1_ini_10km > 36 & s1_ini_10km <= 72] - 36 # Circularity
if(any(terra::values(s1_ini_10km) > 72, na.rm = T)){
  s1_ini_10km[s1_ini_10km > 72] <- s1_ini_10km[s1_ini_10km > 72] - 72
} # Circularity
s1_ini_10km <- terra::subst(s1_ini_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_ini_10km, filename = paste0(out,'/one_s1_ini_10km.tif'), overwrite = T)
rm(s1_ini_10km); gc(TRUE)

# End date
if(!file.exists(paste0(out,'/one_s1_end_raw_10km.tif'))){
  s1_end <- terra::rast(paste0(root,'/phenology/phenoe1_v03.tif'))
  s1_end[s1_end > 250] <- NA
  s1_end[n_seasons == 2] <- NA
  s1_end_10km <- terra::resample(x = s1_end, y = tmp_10km, method = 'near', threads = T)
  terra::writeRaster(x = s1_end_10km, filename = paste0(out,'/one_s1_end_raw_10km.tif'), overwrite = T)
  rm(s1_end); gc(TRUE)
}
s1_end_10km <- terra::rast(paste0(out,'/one_s1_end_raw_10km.tif'))
s1_end_10km[s1_end_10km > 36 & s1_end_10km <= 72] <- s1_end_10km[s1_end_10km > 36 & s1_end_10km <= 72] - 36
if(any(terra::values(s1_end_10km) > 72, na.rm = T)){
  s1_end_10km[s1_end_10km > 72] <- s1_end_10km[s1_end_10km > 72] - 72
}
s1_end_10km <- terra::subst(s1_end_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_end_10km, filename = paste0(out,'/one_s1_end_10km.tif'), overwrite = T)
rm(s1_end_10km); gc(TRUE)

# ------------------------------------------ #
# Two seasons - Season 1 pre-processing
# ------------------------------------------ #
# Start date
if(!file.exists(paste0(out,'/two_s1_ini_raw_10km.tif'))){
  s1_ini <- terra::rast(paste0(root,'/phenology/phenos1_v03.tif'))
  s1_ini[s1_ini > 250] <- NA # Remove missing data
  s1_ini[n_seasons == 1] <- NA
  s1_ini_10km <- terra::resample(x = s1_ini, y = tmp_10km, method = 'near', threads = T)
  rm(s1_ini); gc(TRUE)
  terra::writeRaster(x = s1_ini_10km, filename = paste0(out,'/two_s1_ini_raw_10km.tif'), overwrite = T)
}
s1_ini_10km <- terra::rast(paste0(out,'/two_s1_ini_raw_10km.tif'))
s1_ini_10km[s1_ini_10km > 36 & s1_ini_10km <= 72] <- s1_ini_10km[s1_ini_10km > 36 & s1_ini_10km <= 72] - 36 # Circularity
if(any(terra::values(s1_ini_10km) > 72, na.rm = T)){
  s1_ini_10km[s1_ini_10km > 72] <- s1_ini_10km[s1_ini_10km > 72] - 72
} # Circularity
s1_ini_10km <- terra::subst(s1_ini_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_ini_10km, filename = paste0(out,'/two_s1_ini_10km.tif'), overwrite = T)
rm(s1_ini_10km); gc(TRUE)

# End date
if(!file.exists(paste0(out,'/two_s1_end_raw_10km.tif'))){
  s1_end <- terra::rast(paste0(root,'/phenology/phenoe1_v03.tif'))
  s1_end[s1_end > 250] <- NA
  s1_end[n_seasons == 1] <- NA
  s1_end_10km <- terra::resample(x = s1_end, y = tmp_10km, method = 'near', threads = T)
  terra::writeRaster(x = s1_end_10km, filename = paste0(out,'/two_s1_end_raw_10km.tif'), overwrite = T)
  rm(s1_end); gc(TRUE)
}
s1_end_10km <- terra::rast(paste0(out,'/two_s1_end_raw_10km.tif'))
s1_end_10km[s1_end_10km > 36 & s1_end_10km <= 72] <- s1_end_10km[s1_end_10km > 36 & s1_end_10km <= 72] - 36
if(any(terra::values(s1_end_10km) > 72, na.rm = T)){
  s1_end_10km[s1_end_10km > 72] <- s1_end_10km[s1_end_10km > 72] - 72
}
s1_end_10km <- terra::subst(s1_end_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_end_10km, filename = paste0(out,'/two_s1_end_10km.tif'), overwrite = T)
rm(s1_end_10km); gc(TRUE)

# ------------------------------------------ #
# Two seasons - Season 2 pre-processing
# ------------------------------------------ #
# Start date
if(!file.exists(paste0(out,'/two_s2_ini_raw_10km.tif'))){
  s2_ini <- terra::rast(paste0(root,'/phenology/phenos2_v03.tif'))
  s2_ini[s2_ini > 250] <- NA # Remove missing data
  s2_ini[n_seasons == 1] <- NA
  s2_ini_10km <- terra::resample(x = s2_ini, y = tmp_10km, method = 'near', threads = T)
  rm(s2_ini); gc(TRUE)
  terra::writeRaster(x = s2_ini_10km, filename = paste0(out,'/two_s2_ini_raw_10km.tif'), overwrite = T)
}
s2_ini_10km <- terra::rast(paste0(out,'/two_s2_ini_raw_10km.tif'))
s2_ini_10km[s2_ini_10km > 36 & s2_ini_10km <= 72] <- s2_ini_10km[s2_ini_10km > 36 & s2_ini_10km <= 72] - 36 # Circularity
if(any(terra::values(s2_ini_10km) > 72, na.rm = T)){
  s2_ini_10km[s2_ini_10km > 72] <- s2_ini_10km[s2_ini_10km > 72] - 72
} # Circularity
s2_ini_10km <- terra::subst(s2_ini_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s2_ini_10km, filename = paste0(out,'/two_s2_ini_10km.tif'), overwrite = T)
rm(s2_ini_10km); gc(TRUE)

# End date
if(!file.exists(paste0(out,'/two_s2_end_raw_10km.tif'))){
  s2_end <- terra::rast(paste0(root,'/phenology/phenoe2_v03.tif'))
  s2_end[s2_end > 250] <- NA
  s2_end[n_seasons == 1] <- NA
  s2_end_10km <- terra::resample(x = s2_end, y = tmp_10km, method = 'near', threads = T)
  terra::writeRaster(x = s2_end_10km, filename = paste0(out,'/two_s2_end_raw_10km.tif'), overwrite = T)
  rm(s2_end); gc(TRUE)
}
s2_end_10km <- terra::rast(paste0(out,'/two_s2_end_raw_10km.tif'))
s2_end_10km[s2_end_10km > 36 & s2_end_10km <= 72] <- s2_end_10km[s2_end_10km > 36 & s2_end_10km <= 72] - 36
if(any(terra::values(s2_end_10km) > 72, na.rm = T)){
  s2_end_10km[s2_end_10km > 72] <- s2_end_10km[s2_end_10km > 72] - 72
}
s2_end_10km <- terra::subst(s2_end_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s2_end_10km, filename = paste0(out,'/two_s2_end_10km.tif'), overwrite = T)
rm(s2_end_10km); gc(TRUE)

# ------------------------------------------ #
# Growing season in two years
# ------------------------------------------ #

# One season - Season 1 pre-processing
s1_ini <- terra::rast(paste0(out,'/one_s1_ini_10km.tif'))
s1_end <- terra::rast(paste0(out,'/one_s1_end_10km.tif'))
# 5 km
s1_ini_05km <- terra::resample(x = s1_ini, y = tmp_05km, method = 'near', threads = T)
s1_end_05km <- terra::resample(x = s1_end, y = tmp_05km, method = 'near', threads = T)
cnd <- ((s1_end_05km - s1_ini_05km) <= 0)
s1_ini_05km[cnd] <- NA
s1_end_05km[cnd] <- NA
terra::writeRaster(x = s1_ini_05km, filename = paste0(out,'/one_s1_ini_5km.tif'), overwrite = T, gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
terra::writeRaster(x = s1_end_05km, filename = paste0(out,'/one_s1_end_5km.tif'), overwrite = T, gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
rm(s1_ini_05km, s1_end_05km); gc(T)
# 10 km
s1_ini_10km <- terra::resample(x = s1_ini, y = tmp_10km, method = 'near', threads = T)
s1_end_10km <- terra::resample(x = s1_end, y = tmp_10km, method = 'near', threads = T)
cnd <- ((s1_end_10km - s1_ini_10km) <= 0)
s1_ini_10km[cnd] <- NA
s1_end_10km[cnd] <- NA
terra::writeRaster(x = s1_ini_10km, filename = paste0(out,'/one_s1_ini_10km.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
terra::writeRaster(x = s1_end_10km, filename = paste0(out,'/one_s1_end_10km.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
rm(s1_end, s1_end_05km, s1_end_10km); gc(T)
# Two seasons - Season 1 pre-processing
s1_ini <- terra::rast(paste0(out,'/two_s1_ini.tif'))
s1_ini_05km <- terra::resample(x = s1_ini, y = tmp_05km, method = 'near', threads = T, filename = paste0(out,'/two_s1_ini_5km.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
s1_ini_10km <- terra::resample(x = s1_ini, y = tmp_10km, method = 'near', threads = T, filename = paste0(out,'/two_s1_ini_10km.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
rm(s1_ini, s1_ini_05km, s1_ini_10km); gc(T)
s1_end <- terra::rast(paste0(out,'/two_s1_end.tif'))
s1_end_05km <- terra::resample(x = s1_end, y = tmp_05km, method = 'near', threads = T, filename = paste0(out,'/two_s1_end_5km.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
s1_end_10km <- terra::resample(x = s1_end, y = tmp_10km, method = 'near', threads = T, filename = paste0(out,'/two_s1_end_10km.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
rm(s1_end, s1_end_05km, s1_end_10km); gc(T)
# Two seasons - Season 2 pre-processing
s2_ini <- terra::rast(paste0(out,'/two_s2_ini.tif'))
s2_ini_05km <- terra::resample(x = s2_ini, y = tmp_05km, method = 'near', threads = T, filename = paste0(out,'/two_s2_ini_5km.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
s2_ini_10km <- terra::resample(x = s2_ini, y = tmp_10km, method = 'near', threads = T, filename = paste0(out,'/two_s2_ini_10km.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
rm(s2_ini, s2_ini_05km, s2_ini_10km); gc(T)
s2_end <- terra::rast(paste0(out,'/two_s2_end.tif'))
s2_end_05km <- terra::resample(x = s2_end, y = tmp_05km, method = 'near', threads = T, filename = paste0(out,'/two_s2_end_5km.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
s2_end_10km <- terra::resample(x = s2_end, y = tmp_10km, method = 'near', threads = T, filename = paste0(out,'/two_s2_end_10km.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
rm(s2_end, s2_end_05km, s2_end_10km); gc(T)


## Precipitation template
tmp <- terra::rast(list.files(path = '//CATALOGUE/WFP_ClimateRiskPr1/1.Data/Chirps', pattern = '.tif$', full.names = T)[1])
# tmp <- terra::crop(tmp, terra::ext(colombia))
tmp[tmp == -9999] <- NA

## Re-sampling growing season information to 5 km
s1_ini_5km <- terra::resample(x = s1_ini, y = tmp, method = 'near', threads = T)
s1_end_5km <- terra::resample(x = s1_end, y = tmp, method = 'near', threads = T)

s1_ini_5km[(s1_end_5km - s1_ini_5km) <= 0] <- NA
s1_end_5km[(s1_end_5km - s1_ini_5km) <= 0] <- NA

# if NSeasons == 1
# if Season %in% 1:365 then terra::rapp(...)
# else Season c(Start:365,1:End) then terra::rapp(...)
# sum the two rasters
# if NSeasons == 2
# for i in 1:2
#   if Season %in% 1:365 then terra::rapp(...)
#   else Season c(Start:365,1:End) then terra::rapp(...)
#   sum the two rasters






## Results

options(warn = -1, scipen = 999)
library(pacman)
pacman::p_load(terra, tidyverse)
fls <- list.files(path = 'E:/Harold/Maps', pattern = '.tif$', full.names = T)
p95 <- fls |> purrr::map(.f = function(fl){
  r <- terra::rast(fl)
  idx <- r[[1]]
  whn <- r[[2]]
  return(list(idx, whn))
})
idx <- p95 |> purrr::map(1) |> terra::rast()
whn <- p95 |> purrr::map(2) |> terra::rast()

idx_df <- idx |> terra::as.data.frame(xy = T, cell = T)
whn_df <- whn |> terra::as.data.frame(xy = T, cell = T)

# ----------------------------------------------------------------- #
# Index
# ----------------------------------------------------------------- #

plot(x = 1981:2023,
     y = (as.numeric(idx_df[1,4:ncol(idx_df)])),
     ty = 'l') # ylim = c(0, 202), 
for(i in 2:80){ # nrow(idx_df)
  lines(x = 1981:2023, y = (as.numeric(idx_df[i,4:ncol(idx_df)])), col = 'red')
}
plot(x = 1981:2023,
     y = scale(as.numeric(idx_df[1,4:ncol(idx_df)])),
     ylim = c(-3,3),
     ty = 'l') # ylim = c(0, 202), 
for(i in 2:80){ # nrow(idx_df)
  lines(x = 1981:2023, y = scale(as.numeric(idx_df[i,4:ncol(idx_df)])), col = 'red')
}

# ----------------------------------------------------------------- #
# Date
# ----------------------------------------------------------------- #

plot(x = 1981:2023,
     y = (as.numeric(whn_df[1,4:ncol(whn_df)])),
     ty = 'l') # ylim = c(0, 202), 
for(i in 2:80){ # nrow(idx_df)
  lines(x = 1981:2023, y = (as.numeric(whn_df[i,4:ncol(whn_df)])), col = 'red')
}
plot(x = 1981:2023,
     y = scale(as.numeric(whn_df[1,4:ncol(whn_df)])),
     ylim = c(-3,3),
     ty = 'l') # ylim = c(0, 202), 
for(i in 2:80){ # nrow(idx_df)
  lines(x = 1981:2023, y = scale(as.numeric(whn_df[i,4:ncol(whn_df)])), col = 'red')
}

milk_data %>%
  group_by(Area) %>%
  arrange(Area, Year) %>%
  mutate(Rolling_Median_Production = zoo::rollapply(Value, width = 5, FUN = median, fill = NA, align = "right"),
         Loess_Fitted_Production = loess(Value ~ Year, span = 0.4)$fitted,
         Residuals = Value - Loess_Fitted_Production,
         Lag1_Residuals = lag(Residuals),
         Cooksd = ifelse(!is.na(Lag1_Residuals), cooks.distance(lm(Residuals ~ Lag1_Residuals)), NA),
         Shock_Point = ifelse(Cooksd > 0.1 & Residuals < 0, TRUE, FALSE)))