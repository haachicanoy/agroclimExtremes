# ------------------------------------------ #
# Prepare growing seasons
# By: Harold Achicanoy
# WUR & ABC
# Dec. 2023
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra))

# Directories
root <- 'E:/Harold/Maps' # root <- 'D:/Data/Maps'
out  <- paste0(root,'/results'); dir.create(out, F, T)

# Mapping dekads to day-of-year units
yrs <- 1981:2023
doy <- lapply(X = yrs, FUN = function(yr){
  dks <- seq(from = as.Date(paste0(yr,'-01-01')), to = as.Date(paste0(yr,'-12-31')), length.out = 36)
  doy <- data.frame(dkd = lubridate::yday(dks))
  return(doy)
}) |> dplyr::bind_cols()
doy <- apply(X = doy, MARGIN = 1, FUN = function(x){round(mean(x))})
dks_map <- data.frame(day = doy, dekad = 1:36); rm(yrs, doy)
dks_map$day[dks_map$dekad == 36] <- 1

# Load phenology data
n_seasons <- terra::rast(paste0(root,'/phenology/phenonseasons_v03.tif')); gc(TRUE)

# ------------------------------------------ #
# One season - Season 1 pre-processing
# ------------------------------------------ #
# Start date
s1_ini <- terra::rast(paste0(root,'/phenology/phenos1_v03.tif'))
s1_ini[s1_ini > 250] <- NA # Remove missing data
s1_ini[n_seasons != 1] <- NA
s1_ini[s1_ini > 36 & s1_ini <= 72] <- s1_ini[s1_ini > 36 & s1_ini <= 72] - 36 # Circularity
s1_ini[s1_ini > 72] <- s1_ini[s1_ini > 72] - 72 # Circularity
for(i in 1:nrow(dks_map)){ # Dekad to day-of-year
  s1_ini[s1_ini == i] <- dks_map$day[i]
  gc(TRUE)
}; rm(i)
terra::writeRaster(s1_ini, filename = paste0(out,'/one_s1_ini.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
rm(s1_ini); gc(TRUE)

# End date
s1_end <- terra::rast(paste0(root,'/phenology/phenoe1_v03.tif'))
s1_end[s1_end > 250] <- NA
s1_end[n_seasons != 1] <- NA
s1_end[s1_end > 36 & s1_end <= 72] <- s1_end[s1_end > 36 & s1_end <= 72] - 36
s1_end[s1_end > 72] <- s1_end[s1_end > 72] - 72
for(i in 1:nrow(dks_map)){
  s1_end[s1_end == i] <- dks_map$day[i]
  gc(TRUE)
}; rm(i)
terra::writeRaster(s1_end, filename = paste0(out,'/one_s1_end.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
rm(s1_end); gc(TRUE)

# ------------------------------------------ #
# Two seasons - Season 1 pre-processing
# ------------------------------------------ #
# Start date
s1_ini <- terra::rast(paste0(root,'/phenology/phenos1_v03.tif'))
s1_ini[s1_ini > 250] <- NA
s1_ini[n_seasons != 2] <- NA
s1_ini[s1_ini > 36 & s1_ini <= 72] <- s1_ini[s1_ini > 36 & s1_ini <= 72] - 36
s1_ini[s1_ini > 72] <- s1_ini[s1_ini > 72] - 72
for(i in 1:nrow(dks_map)){
  s1_ini[s1_ini == i] <- dks_map$day[i]
  gc(TRUE)
}; rm(i)
terra::writeRaster(s1_ini, filename = paste0(out,'/two_s1_ini.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
rm(s1_ini); gc(TRUE)

# End date
s1_end <- terra::rast(paste0(root,'/phenology/phenoe1_v03.tif'))
s1_end[s1_end > 250] <- NA
s1_end[n_seasons != 2] <- NA
s1_end[s1_end > 36 & s1_end <= 72] <- s1_end[s1_end > 36 & s1_end <= 72] - 36
s1_end[s1_end > 72] <- s1_end[s1_end > 72] - 72
for(i in 1:nrow(dks_map)){
  s1_end[s1_end == i] <- dks_map$day[i]
  gc(TRUE)
}; rm(i)
terra::writeRaster(s1_end, filename = paste0(out,'/two_s1_end.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
rm(s1_end); gc(TRUE)

# ------------------------------------------ #
# Two seasons - Season 2 pre-processing
# ------------------------------------------ #
# Start date
s2_ini <- terra::rast(paste0(root,'/phenology/phenos2_v03.tif'))
s2_ini[s2_ini > 250] <- NA
s2_ini[n_seasons != 2] <- NA
s2_ini[s2_ini > 36 & s2_ini <= 72] <- s2_ini[s2_ini > 36 & s2_ini <= 72] - 36
s2_ini[s2_ini > 72] <- s2_ini[s2_ini > 72] - 72
for(i in 1:nrow(dks_map)){
  s2_ini[s2_ini == i] <- dks_map$day[i]
  gc(TRUE)
}; rm(i)
terra::writeRaster(s2_ini, filename = paste0(out,'/two_s2_ini.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
rm(s2_ini); gc(TRUE)

# End date
s2_end <- terra::rast(paste0(root,'/phenology/phenoe2_v03.tif'))
s2_end[s2_end > 250] <- NA
s2_end[n_seasons != 2] <- NA
s2_end[s2_end > 36 & s2_end <= 72] <- s2_end[s2_end > 36 & s2_end <= 72] - 36
s2_end[s2_end > 72] <- s2_end[s2_end > 72] - 72
for(i in 1:nrow(dks_map)){
  s2_end[s2_end == i] <- dks_map$day[i]
  gc(TRUE)
}; rm(i)
terra::writeRaster(s2_end, filename = paste0(out,'/two_s2_end.tif'), gdal = c('COMPRESS=NONE', 'TFW=YES'), datatype = 'INT1U')
rm(s2_end); gc(TRUE)





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