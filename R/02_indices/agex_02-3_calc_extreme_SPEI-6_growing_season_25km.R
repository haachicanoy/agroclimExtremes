# ------------------------------------------ #
# Get extreme SPEI value within growing season
# By: Harold Achicanoy
# WUR & ABC
# Dec. 2023
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,lubridate))

# Root directory
root <- '//CATALOGUE/AgroclimExtremes'

# Source agro-climatic indices functions
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')
# Find minimum
fmin <- function(x){
  x[is.infinite(x)] <- NA
  min(x)
}

# Spatial resolution
res <- 25

# SPEI files
fls <- list.files(path = paste0(root,'/agex_raw_data/monthly_spei_',res,'km'), pattern = '.tif$', full.names = T)
fls <- fls[-(1:5)]

# ------------------------------------------ #
# One season
# ------------------------------------------ #

# There are two possibilities. One-year growing season which indicates that
# start and ending dates are within the same year. Two-years growing season
# indicating that the start and ending dates cover two consecutive years.

# >>> Processing of One-year growing season
s1_ini <- terra::rast(paste0(root,'/agex_raw_data/agex_one_s1_ini_',res,'km_croplands.tif'))
s1_end <- terra::rast(paste0(root,'/agex_raw_data/agex_one_s1_end_',res,'km_croplands.tif'))
s1_dff <- s1_end - s1_ini
# Selecting the pixels with One-year growing season
s1_dff[s1_end < s1_ini] <- NA
s1_ini[is.na(s1_dff)] <- NA
s1_end[is.na(s1_dff)] <- NA
# Day of the year to month
s1_ini <- terra::app(x = s1_ini, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})
s1_end <- terra::app(x = s1_end, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})

# >>> Indices computation
# Maximum SPEI within growing season
yrs <- 1979:2023
lapply(X = yrs, FUN = function(yr){
  fls_flt <- fls[grep(pattern = paste0('_',yr,'-'), x = fls)] # Files filtered per year
  Spei <- terra::rast(fls_flt); gc(T)
  mxSpei <- terra::rapp(x = Spei, first = s1_ini, last = s1_end, fun = fmin); gc(T)
  mxSpei <- -1 * mxSpei
  names(mxSpei) <- paste0('spei-6_',yr)
  outfile <- paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_intermediate_',res,'km/one_s1_y1_spei-6_',yr,'.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = mxSpei, filename = outfile, overwrite = T); gc(T)
  return('Done.\n')
})

# >>> Processing of Two-years growing season
s1_ini <- terra::rast(paste0(root,'/agex_raw_data/agex_one_s1_ini_',res,'km_croplands.tif'))
s1_end <- terra::rast(paste0(root,'/agex_raw_data/agex_one_s1_end_',res,'km_croplands.tif'))
s1_dff <- s1_end - s1_ini
# Selecting the pixels with Two-years growing season
s1_dff[s1_end > s1_ini] <- NA
s1_ini[is.na(s1_dff)] <- NA
s1_end[is.na(s1_dff)] <- NA
# Day of the year to month
s1_ini <- terra::app(x = s1_ini, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})
s1_end <- terra::app(x = s1_end, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})

# >>> Indices computation
# Maximum SPEI within growing season
yrs <- 1979:2023
lapply(X = yrs, FUN = function(yr){
  if(yr < 2023){
    fls_flt <- fls[unlist(grep2(pattern = c(paste0('_',yr,'-'),paste0('_',yr+1,'-')), x = fls))] # To implement
    Spei <- terra::rast(fls_flt)
    mxSpei <- terra::rapp(x = Spei, first = s1_ini, last = s1_end+12, fun = fmin)
    mxSpei <- -1 * mxSpei
    names(mxSpei) <- paste0('spei-6_',yr)
    outfile <- paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_intermediate_',res,'km/one_s1_y2_spei-6_',yr,'.tif')
    dir.create(dirname(outfile),F,T)
    terra::writeRaster(x = mxSpei, filename = outfile, overwrite = T); gc(T)
    return('Done.\n')
  } else {
    break
    return('Done.\n')
  }
})

# >>> Merging results of One-year and Two-years growing seasons
# Maximum number of consecutive dry days
yrs <- 1980:2022
lapply(X = yrs, FUN = function(yr){
  fls <- list.files(pattern = paste0('_',yr,'.tif$'), path = paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_intermediate_',res,'km'), full.names = T)
  mxSpei <- lapply(fls, terra::rast)
  mxSpei <- terra::merge(x = mxSpei[[1]], y = mxSpei[[2]])
  outfile <- paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_',res,'km/one_s1_spei-6_',yr,'.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = mxSpei, filename = outfile, overwrite = T)
})

# ------------------------------------------ #
# Two seasons
# ------------------------------------------ #

# There are two possibilities. One-year growing season which indicates that
# start and ending dates are within the same year. Two-years growing season
# indicating that the start and ending dates cover two consecutive years.
# But in this case, it applies for both growing seasons.

## Season 1
# >>> Processing of One-year growing season
s1_ini <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s1_ini_',res,'km_croplands.tif'))
s1_end <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s1_end_',res,'km_croplands.tif'))
s1_dff <- s1_end - s1_ini
# Selecting the pixels with One-year growing season
s1_dff[s1_end < s1_ini] <- NA
s1_ini[is.na(s1_dff)] <- NA
s1_end[is.na(s1_dff)] <- NA
# Day of the year to month
s1_ini <- terra::app(x = s1_ini, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})
s1_end <- terra::app(x = s1_end, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})

# >>> Indices computation
# Maximum SPEI within growing season
yrs <- 1979:2023
lapply(X = yrs, FUN = function(yr){
  fls_flt <- fls[grep(pattern = paste0('_',yr,'-'), x = fls)] # Files filtered per year
  Spei <- terra::rast(fls_flt); gc(T)
  mxSpei <- terra::rapp(x = Spei, first = s1_ini, last = s1_end, fun = fmin); gc(T)
  mxSpei <- -1 * mxSpei
  names(mxSpei) <- paste0('spei-6_',yr)
  outfile <- paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_intermediate_',res,'km/two_s1_y1_spei-6_',yr,'.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = mxSpei, filename = outfile, overwrite = T); gc(T)
  return('Done.\n')
})

# >>> Processing of Two-years growing season
s1_ini <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s1_ini_',res,'km_croplands.tif'))
s1_end <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s1_end_',res,'km_croplands.tif'))
s1_dff <- s1_end - s1_ini
# Selecting the pixels with Two-years growing season
s1_dff[s1_end > s1_ini] <- NA
s1_ini[is.na(s1_dff)] <- NA
s1_end[is.na(s1_dff)] <- NA
# Day of the year to month
s1_ini <- terra::app(x = s1_ini, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})
s1_end <- terra::app(x = s1_end, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})

# >>> Indices computation
# Maximum SPEI within growing season
yrs <- 1979:2023
lapply(X = yrs, FUN = function(yr){
  if(yr < 2023){
    fls_flt <- fls[unlist(grep2(pattern = c(paste0('_',yr,'-'),paste0('_',yr+1,'-')), x = fls))]
    Spei <- terra::rast(fls_flt)
    mxSpei <- terra::rapp(x = Spei, first = s1_ini, last = s1_end+12, fun = fmin)
    mxSpei <- -1 * mxSpei
    names(mxSpei) <- paste0('spei-6_',yr)
    outfile <- paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_intermediate_',res,'km/two_s1_y2_spei-6_',yr,'.tif')
    dir.create(dirname(outfile),F,T)
    terra::writeRaster(x = mxSpei, filename = outfile, overwrite = T); gc(T)
    return('Done.\n')
  } else {
    break
    return('Done.\n')
  }
})

# >>> Merging results of One-year and Two-years growing seasons
# Maximum SPEI within growing season
yrs <- 1980:2022
lapply(X = yrs, FUN = function(yr){
  fls <- list.files(pattern = paste0('two_s1_y[0-9]_spei-6_',yr,'.tif'), path = paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_intermediate_',res,'km'), full.names = T)
  mxSpei <- lapply(fls, terra::rast)
  mxSpei <- terra::merge(x = mxSpei[[1]], y = mxSpei[[2]])
  outfile <- paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_',res,'km/two_s1_spei-6_',yr,'.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = mxSpei, filename = outfile, overwrite = T)
})

## Season 2
# >>> Processing of One-year growing season
s2_ini <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s2_ini_',res,'km_croplands.tif'))
s2_end <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s2_end_',res,'km_croplands.tif'))
s2_dff <- s2_end - s2_ini
# Selecting the pixels with One-year growing season
s2_dff[s2_end < s2_ini] <- NA
s2_ini[is.na(s2_dff)] <- NA
s2_end[is.na(s2_dff)] <- NA
# Day of the year to month
s2_ini <- terra::app(x = s2_ini, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})
s2_end <- terra::app(x = s2_end, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})

# >>> Indices computation
# Maximum SPEI within growing season
yrs <- 1979:2023
lapply(X = yrs, FUN = function(yr){
  fls_flt <- fls[grep(pattern = paste0('_',yr,'-'), x = fls)] # Files filtered per year
  Spei <- terra::rast(fls_flt); gc(T)
  mxSpei <- terra::rapp(x = Spei, first = s2_ini, last = s2_end, fun = fmin); gc(T)
  mxSpei <- -1 * mxSpei
  names(mxSpei) <- paste0('spei-6_',yr)
  outfile <- paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_intermediate_',res,'km/two_s2_y1_spei-6_',yr,'.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = mxSpei, filename = outfile, overwrite = T); gc(T)
  return('Done.\n')
})

# >>> Processing of Two-years growing season
s2_ini <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s2_ini_',res,'km_croplands.tif'))
s2_end <- terra::rast(paste0(root,'/agex_raw_data/agex_two_s2_end_',res,'km_croplands.tif'))
s2_dff <- s2_end - s2_ini
# Selecting the pixels with Two-years growing season
s2_dff[s2_end > s2_ini] <- NA
s2_ini[is.na(s2_dff)] <- NA
s2_end[is.na(s2_dff)] <- NA
# Day of the year to month
s2_ini <- terra::app(x = s2_ini, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})
s2_end <- terra::app(x = s2_end, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})

# >>> Indices computation
# Maximum SPEI within growing season
yrs <- 1979:2023
lapply(X = yrs, FUN = function(yr){
  if(yr < 2023){
    fls_flt <- fls[unlist(grep2(pattern = c(paste0('_',yr,'-'),paste0('_',yr+1,'-')), x = fls))]
    Spei <- terra::rast(fls_flt)
    mxSpei <- terra::rapp(x = Spei, first = s2_ini, last = s2_end+12, fun = fmin)
    mxSpei <- -1 * mxSpei
    names(mxSpei) <- paste0('spei-6_',yr)
    outfile <- paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_intermediate_',res,'km/two_s2_y2_spei-6_',yr,'.tif')
    dir.create(dirname(outfile),F,T)
    terra::writeRaster(x = mxSpei, filename = outfile, overwrite = T); gc(T)
    return('Done.\n')
  } else {
    break
    return('Done.\n')
  }
})

# >>> Merging results of One-year and Two-years growing seasons
# Maximum SPEI within growing season
yrs <- 1980:2022
lapply(X = yrs, FUN = function(yr){
  fls <- list.files(pattern = paste0('two_s2_y[0-9]_spei-6_',yr,'.tif'), path = paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_intermediate_',res,'km'), full.names = T)
  mxSpei <- lapply(fls, terra::rast)
  mxSpei <- terra::merge(x = mxSpei[[1]], y = mxSpei[[2]])
  outfile <- paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_',res,'km/two_s2_spei-6_',yr,'.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = mxSpei, filename = outfile, overwrite = T)
})

# # Creating raster stack files
# index  <- 'spei-6'
# gs     <- 'two'
# season <- 2
# res    <- 25
# pth    <- paste0(root,'/agex_indices/agex_',index,'/agex_',index,'_',res,'km')
# 
# fls <- list.files(path = pth,
#                   pattern = paste0(gs,'_s',season,'_',index,'_*.*.tif$'),
#                   full.names = T)
# r <- terra::rast(fls)
# terra::writeRaster(x = r, filename = paste0(pth,'/',gs,'_s',season,'_',index,'_25km.tif'), overwrite = T)
