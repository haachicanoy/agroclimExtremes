# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Calculate monthly SPEI-6 at 0.1°
# By: Harold Achicanoy, Cesar Saavedra
# WUR & ABC
# Created in February 2024
# Modified in January 2026
# --------------------------------------------------------------- #

# R options and user-defined functions
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,SPEI))

# Root directory
ddir <- '//CATALOGUE/WFP_ClimateRiskPr1'
root <- '//CATALOGUE/AgroclimExtremes'

# List water balance files
blc_fls <- list.files(path = paste0(root,'/agex_raw_data/monthly_balance'), pattern = '.tif$', full.names = T)
nms_fls <- gsub('.tif','', gsub('bal-', '', basename(blc_fls)))

# Out files
out_pth <- paste0(root,'/agex_raw_data/monthly_spei')
out_fls <- paste0(out_pth,'/spei-',scl,'_',nms_fls,'.tif')

# Load water balance files
blc <- terra::rast(blc_fls); names(blc) <- nms_fls; gc(T)

calc_spei <- function(x){
  x <- ts(data = x, start = c(1979,1), frequency = 12)
  Spei <- SPEI::spei(data = x, scale = 6, na.rm = T)$fitted |> as.numeric()
  return(Spei)
}

# Calculate SPEI
# Spei <- terra::app(x = blc, fun = calc_spei, cores = 20)
Spei <- terra::app(x = blc, fun = function(i, ff) ff(i), cores = 20, ff = calc_spei)
names(Spei) <- nms_fls
terra::writeRaster(Spei, out_fls, overwrite = T)
