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

# Template rasters
tmp_10km <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5.tif')

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
n_seasons <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_phenology/phenonseasons_v03.tif')); gc(T)
n_seasons <- terra::resample(x = n_seasons, y = tmp_10km, method = 'near', threads = T); gc(T) # Growing seasons per site at 10 km resolution

# ------------------------------------------ #
# One season - Season 1 pre-processing
# ------------------------------------------ #
# Start date
if(!file.exists(paste0(aux,'/agex_one_s1_ini_raw_10km.tif'))){
  s1_ini <- terra::rast(paste0(aux,'/phenos1_v03.tif'))
  s1_ini[s1_ini > 250] <- NA # Remove missing data
  s1_ini_10km <- terra::resample(x = s1_ini, y = tmp_10km, method = 'near', threads = T)
  s1_ini_10km[n_seasons == 2] <- NA
  terra::writeRaster(x = s1_ini_10km, filename = paste0(aux,'/agex_one_s1_ini_raw_10km.tif'), overwrite = T)
  rm(s1_ini); gc(T)
}
s1_ini_10km <- terra::rast(paste0(aux,'/agex_one_s1_ini_raw_10km.tif'))
s1_ini_10km[s1_ini_10km > 36 & s1_ini_10km <= 72] <- s1_ini_10km[s1_ini_10km > 36 & s1_ini_10km <= 72] - 36 # Circularity
if(any(terra::values(s1_ini_10km) > 72, na.rm = T)){
  s1_ini_10km[s1_ini_10km > 72] <- s1_ini_10km[s1_ini_10km > 72] - 72
} # Circularity
s1_ini_10km <- terra::subst(s1_ini_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_ini_10km, filename = paste0(out,'/agex_one_s1_ini_10km.tif'), overwrite = T)
rm(s1_ini_10km); gc(T)

# End date
if(!file.exists(paste0(aux,'/agex_one_s1_end_raw_10km.tif'))){
  s1_end <- terra::rast(paste0(aux,'/phenoe1_v03.tif'))
  s1_end[s1_end > 250] <- NA
  s1_end_10km <- terra::resample(x = s1_end, y = tmp_10km, method = 'near', threads = T)
  s1_end_10km[n_seasons == 2] <- NA
  terra::writeRaster(x = s1_end_10km, filename = paste0(aux,'/agex_one_s1_end_raw_10km.tif'), overwrite = T)
  rm(s1_end); gc(T)
}
s1_end_10km <- terra::rast(paste0(aux,'/agex_one_s1_end_raw_10km.tif'))
s1_end_10km[s1_end_10km > 36 & s1_end_10km <= 72] <- s1_end_10km[s1_end_10km > 36 & s1_end_10km <= 72] - 36
if(any(terra::values(s1_end_10km) > 72, na.rm = T)){
  s1_end_10km[s1_end_10km > 72] <- s1_end_10km[s1_end_10km > 72] - 72
}
s1_end_10km <- terra::subst(s1_end_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_end_10km, filename = paste0(out,'/agex_one_s1_end_10km.tif'), overwrite = T)
rm(s1_end_10km); gc(T)

# ------------------------------------------ #
# Two seasons - Season 1 pre-processing
# ------------------------------------------ #
# Start date
if(!file.exists(paste0(aux,'/agex_two_s1_ini_raw_10km.tif'))){
  s1_ini <- terra::rast(paste0(aux,'/phenos1_v03.tif'))
  s1_ini[s1_ini > 250] <- NA # Remove missing data
  s1_ini_10km <- terra::resample(x = s1_ini, y = tmp_10km, method = 'near', threads = T)
  s1_ini_10km[n_seasons == 1] <- NA
  terra::writeRaster(x = s1_ini_10km, filename = paste0(aux,'/agex_two_s1_ini_raw_10km.tif'), overwrite = T)
  rm(s1_ini); gc(T)
}
s1_ini_10km <- terra::rast(paste0(aux,'/agex_two_s1_ini_raw_10km.tif'))
s1_ini_10km[s1_ini_10km > 36 & s1_ini_10km <= 72] <- s1_ini_10km[s1_ini_10km > 36 & s1_ini_10km <= 72] - 36 # Circularity
if(any(terra::values(s1_ini_10km) > 72, na.rm = T)){
  s1_ini_10km[s1_ini_10km > 72] <- s1_ini_10km[s1_ini_10km > 72] - 72
} # Circularity
s1_ini_10km <- terra::subst(s1_ini_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_ini_10km, filename = paste0(out,'/agex_two_s1_ini_10km.tif'), overwrite = T)
rm(s1_ini_10km); gc(T)

# End date
if(!file.exists(paste0(aux,'/agex_two_s1_end_raw_10km.tif'))){
  s1_end <- terra::rast(paste0(aux,'/phenoe1_v03.tif'))
  s1_end[s1_end > 250] <- NA
  s1_end_10km <- terra::resample(x = s1_end, y = tmp_10km, method = 'near', threads = T)
  s1_end_10km[n_seasons == 1] <- NA
  terra::writeRaster(x = s1_end_10km, filename = paste0(aux,'/agex_two_s1_end_raw_10km.tif'), overwrite = T)
  rm(s1_end); gc(T)
}
s1_end_10km <- terra::rast(paste0(aux,'/agex_two_s1_end_raw_10km.tif'))
s1_end_10km[s1_end_10km > 36 & s1_end_10km <= 72] <- s1_end_10km[s1_end_10km > 36 & s1_end_10km <= 72] - 36
if(any(terra::values(s1_end_10km) > 72, na.rm = T)){
  s1_end_10km[s1_end_10km > 72] <- s1_end_10km[s1_end_10km > 72] - 72
}
s1_end_10km <- terra::subst(s1_end_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s1_end_10km, filename = paste0(out,'/agex_two_s1_end_10km.tif'), overwrite = T)
rm(s1_end_10km); gc(T)

# ------------------------------------------ #
# Two seasons - Season 2 pre-processing
# ------------------------------------------ #
# Start date
if(!file.exists(paste0(aux,'/agex_two_s2_ini_raw_10km.tif'))){
  s2_ini <- terra::rast(paste0(aux,'/phenos2_v03.tif'))
  s2_ini[s2_ini > 250] <- NA # Remove missing data
  s2_ini_10km <- terra::resample(x = s2_ini, y = tmp_10km, method = 'near', threads = T)
  s2_ini_10km[n_seasons == 1] <- NA
  terra::writeRaster(x = s2_ini_10km, filename = paste0(aux,'/agex_two_s2_ini_raw_10km.tif'), overwrite = T)
  rm(s2_ini); gc(T)
}
s2_ini_10km <- terra::rast(paste0(aux,'/agex_two_s2_ini_raw_10km.tif'))
s2_ini_10km[s2_ini_10km > 36 & s2_ini_10km <= 72] <- s2_ini_10km[s2_ini_10km > 36 & s2_ini_10km <= 72] - 36 # Circularity
if(any(terra::values(s2_ini_10km) > 72, na.rm = T)){
  s2_ini_10km[s2_ini_10km > 72] <- s2_ini_10km[s2_ini_10km > 72] - 72
} # Circularity
s2_ini_10km <- terra::subst(s2_ini_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s2_ini_10km, filename = paste0(out,'/agex_two_s2_ini_10km.tif'), overwrite = T)
rm(s2_ini_10km); gc(T)

# End date
if(!file.exists(paste0(aux,'/agex_two_s2_end_raw_10km.tif'))){
  s2_end <- terra::rast(paste0(aux,'/phenoe2_v03.tif'))
  s2_end[s2_end > 250] <- NA
  s2_end_10km <- terra::resample(x = s2_end, y = tmp_10km, method = 'near', threads = T)
  s2_end_10km[n_seasons == 1] <- NA
  terra::writeRaster(x = s2_end_10km, filename = paste0(aux,'/agex_two_s2_end_raw_10km.tif'), overwrite = T)
  rm(s2_end); gc(T)
}
s2_end_10km <- terra::rast(paste0(aux,'/agex_two_s2_end_raw_10km.tif'))
s2_end_10km[s2_end_10km > 36 & s2_end_10km <= 72] <- s2_end_10km[s2_end_10km > 36 & s2_end_10km <= 72] - 36
if(any(terra::values(s2_end_10km) > 72, na.rm = T)){
  s2_end_10km[s2_end_10km > 72] <- s2_end_10km[s2_end_10km > 72] - 72
}
s2_end_10km <- terra::subst(s2_end_10km, from = dks_map$dekad, to = dks_map$day)
terra::writeRaster(s2_end_10km, filename = paste0(out,'/agex_two_s2_end_10km.tif'), overwrite = T)
rm(s2_end_10km); gc(T)

## ------------------------------------------ ##
## Filtering of growing seasons within croplands
## ------------------------------------------ ##

# One season
s1_ini_10km <- terra::rast(paste0(out,'/agex_one_s1_ini_10km.tif'))
croplnds    <- terra::rast(paste0(out,'/agex_croplands_mapspam_10km.tif'))
s1_ini_10km_msk <- terra::mask(x = s1_ini_10km, mask = croplnds)
terra::writeRaster(x = s1_ini_10km_msk, filename = paste0(out,'/agex_one_s1_ini_10km_croplands.tif'))
s1_end_10km <- terra::rast(paste0(out,'/agex_one_s1_end_10km.tif'))
s1_end_10km_msk <- terra::mask(x = s1_end_10km, mask = croplnds)
terra::writeRaster(x = s1_end_10km_msk, filename = paste0(out,'/agex_one_s1_end_10km_croplands.tif'))
rm(s1_ini_10km, s1_ini_10km_msk, s1_end_10km, s1_end_10km_msk); gc(T)

# Two seasons - Season 1
s1_ini_10km <- terra::rast(paste0(out,'/agex_two_s1_ini_10km.tif'))
s1_ini_10km_msk <- terra::mask(x = s1_ini_10km, mask = croplnds)
terra::writeRaster(x = s1_ini_10km_msk, filename = paste0(out,'/agex_two_s1_ini_10km_croplands.tif'))
s1_end_10km <- terra::rast(paste0(out,'/agex_two_s1_end_10km.tif'))
s1_end_10km_msk <- terra::mask(x = s1_end_10km, mask = croplnds)
terra::writeRaster(x = s1_end_10km_msk, filename = paste0(out,'/agex_two_s1_end_10km_croplands.tif'))
rm(s1_ini_10km, s1_ini_10km_msk, s1_end_10km, s1_end_10km_msk); gc(T)
# Two seasons - Season 2
s2_ini_10km <- terra::rast(paste0(out,'/agex_two_s2_ini_10km.tif'))
s2_ini_10km_msk <- terra::mask(x = s2_ini_10km, mask = croplnds)
terra::writeRaster(x = s2_ini_10km_msk, filename = paste0(out,'/agex_two_s2_ini_10km_croplands.tif'))
s2_end_10km <- terra::rast(paste0(out,'/agex_two_s2_end_10km.tif'))
s2_end_10km_msk <- terra::mask(x = s2_end_10km, mask = croplnds)
terra::writeRaster(x = s2_end_10km_msk, filename = paste0(out,'/agex_two_s2_end_10km_croplands.tif'))
rm(s2_ini_10km, s2_ini_10km_msk, s2_end_10km, s2_end_10km_msk); gc(T)

# ------------------------------------------ #
# Growing season in two years
# ------------------------------------------ #

# One season - Season 1
s1_ini <- terra::rast(paste0(out,'/agex_one_s1_ini_10km_croplands.tif'))
s1_end <- terra::rast(paste0(out,'/agex_one_s1_end_10km_croplands.tif'))

# Length of the growing season
s1_dff <- s1_end - s1_ini

plot(s1_dff > 0) # growing season within the same year
plot(s1_dff < 0) # growing season in two years

msk_two_yrs <- s1_dff < 0; msk_two_yrs <- msk_two_yrs * 1; msk_two_yrs[msk_two_yrs == 0] <- NA
tst <- terra::mask(x = s1_ini, mask = msk_two_yrs)


vls <- terra::values(s1_dff) |> as.numeric()
vls <- vls[!is.na(vls)]

vls_pos <- vls[vls > 0]
vls_neg <- vls[vls < 0]
summary(vls_pos); summary(vls_neg)

qtl_pos <- stats::quantile(x = vls_pos, probs = 0.95) # c(0.94,0.95,0.96, 0.97,0.98)
qtl_neg <- stats::quantile(x = vls_neg, probs = 0.05) # c(0.02,0.03,0.04, 0.05,0.06)

hist(vls_pos)
abline(v = qtl_pos, col = 2)

hist(vls_neg)
abline(v = qtl_neg, col = 2)

plot(s1_dff > qtl_pos)
plot(s1_dff < qtl_neg)

s1 <- c(s1_ini, s1_end, s1_dff)
s1_dfm <- terra::as.data.frame(x = s1, xy = T, cell = T, na.rm = T)
names(s1_dfm) <- c('cell','x','y','start','end','difference')

rbind(
unique(s1_dfm[s1_dfm$difference == qtl_pos,c('start','end','difference')]),
unique(s1_dfm[s1_dfm$difference > qtl_pos,c('start','end','difference')])
)

plot(s1_dff < 0)
plot(s1_end[s1_dff < 0])

# days in two years
dys <- c(292:365, 1:32)
dys <- c(292:365, 365+1:32)
# Two consecutive years
yrs <- 1979:1980
# Consecutive dates
dts <- seq(from = as.Date(paste0(yrs[1],'-01-01')), to = as.Date(paste0(yrs[2],'-12-31')), by = 'day')


s1_dff[s1_end > s1_ini] <- NA
if(any(terra::values(s1_dff, na.rm = T) > -40)){
  s1_dff[s1_dff > -40] <- NA
}
plot(!is.na(s1_dff))

# Two seasons - Season 1 pre-processing
s1_ini <- terra::rast(paste0(out,'/two_s1_ini_10km.tif'))
s1_end <- terra::rast(paste0(out,'/two_s1_end_10km.tif'))

s1_dff <- s1_end - s1_ini
s1_dff[s1_end > s1_ini] <- NA
if(any(terra::values(s1_dff, na.rm = T) > -40)){
  s1_dff[s1_dff > -40] <- NA
}
plot(!is.na(s1_dff))

# Two seasons - Season 1 pre-processing
s2_ini <- terra::rast(paste0(out,'/two_s2_ini_10km.tif'))
s2_end <- terra::rast(paste0(out,'/two_s2_end_10km.tif'))

s2_dff <- s2_end - s2_ini
s2_dff[s2_end > s2_ini] <- NA
if(any(terra::values(s2_dff, na.rm = T) > -40)){
  s2_dff[s2_dff > -40] <- NA
}
plot(!is.na(s2_dff))

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