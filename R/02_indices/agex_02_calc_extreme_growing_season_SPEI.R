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
root <- '//CATALOGUE/WFP_ClimateRiskPr1'

# Source agro-climatic indices functions
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')

# SPEI files
fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_raw_data/monthly_spei'), pattern = '.tif$', full.names = T)
fls <- fls[-(1:5)]

# ------------------------------------------ #
# One season
# ------------------------------------------ #

# There are two possibilities. One-year growing season which indicates that
# start and ending dates are within the same year. Two-years growing season
# indicating that the start and ending dates cover two consecutive years.

# >>> Processing of One-year growing season
s1_ini <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_ini_10km_croplands.tif'))
s1_end <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_end_10km_croplands.tif'))
s1_dff <- s1_end - s1_ini
# Selecting the pixels with One-year growing season
s1_dff[s1_end < s1_ini] <- NA
s1_ini[is.na(s1_dff)] <- NA
s1_end[is.na(s1_dff)] <- NA
# Day of the year to month
s1_ini <- terra::app(x = s1_ini, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})
s1_end <- terra::app(x = s1_end, function(x){x <- as.Date(x, origin = '2014-01-01'); x <- lubridate::month(x); return(x)})

# Find minimum
fmin <- function(x){
  x[is.infinite(x)] <- NA
  min(x)
}

# >>> Indices computation
# Maximum SPEI within growing season
lapply(X = yrs, FUN = function(yr){
  fls_flt <- fls[grep(pattern = paste0('_',yr,'-'), x = fls)] # Files filtered per year
  Spei <- terra::rast(fls_flt); gc(T)
  mxSpei <- terra::rapp(x = Spei, first = s1_ini, last = s1_end, fun = fmin); gc(T)
  mxSpei <- -1 * mxSpei
  names(mxSpei) <- paste0('spei_',yr)
  terra::writeRaster(x = rsl, filename = paste0(root,'/agroclimExtremes/agex_indices/agex_spei/agex_spei_intermediate/one_s1_y1_spei_',yr,'.tif'), overwrite = T); gc(T)
  return('Done.\n')
})

# >>> Processing of Two-years growing season
s1_ini <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_ini_10km_croplands.tif'))
s1_end <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_one_s1_end_10km_croplands.tif'))
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
lapply(X = yrs, FUN = function(yr){
  if(yr < 2023){
    fls_flt <- prc[unlist(grep2(pattern = c(paste0('_',yr,'-'),paste0('_',yr+1,'-')), x = fls))] # To implement
    Spei <- terra::rast(fls_flt)
    dys <- ifelse(lubridate::leap_year(yr),366,365)
    mxSpei <- terra::rapp(x = Spei, first = s1_ini, last = s1_end+dys, fun = fmin) # to check
    mxSpei <- -1 * mxSpei
    names(mxSpei) <- paste0('spei_',yr)
    terra::writeRaster(x = rsl, filename = paste0(root,'/agroclimExtremes/agex_indices/agex_cdd/agex_cdd_intermediate/one_s1_y2_cdd_',yr,'.tif'), overwrite = T); gc(T)
    return('Done.\n')
  } else {
    break
    return('Done.\n')
  }
})
