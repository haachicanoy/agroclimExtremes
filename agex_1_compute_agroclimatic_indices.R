# ------------------------------------------ #
# Compute agro-climatic indices
# By: Harold Achicanoy
# WUR & ABC
# Dec. 2023
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra))

# Source agro-climatic indices functions
source('.../agex_0_indices.R')

# Precipitation files
prc <- list.files(path = '//CATALOGUE/WFP_ClimateRiskPr1/1.Data/Chirps', pattern = '.tif$', full.names = T)

# ------------------------------------------ #
# One season
# ------------------------------------------ #
# One-year
s1_ini_5km <- terra::rast(paste0(dir,'/one_s1_ini_5km.tif'))
s1_end_5km <- terra::rast(paste0(dir,'/one_s1_end_5km.tif'))

yrs <- 1981:2022
p95_index <- lapply(X = yrs, FUN = function(yr){
  # Precipitation files filtered
  prc_flt <- prc[grep(pattern = yr, x = prc)]
  rnf <- terra::rast(prc_flt)
  rsl <- terra::rapp(x = rnf, first = s1_ini_5km, last = s1_end_5km, fun = p95)
  terra::writeRaster(x = rsl, filename = paste0(root,'/yr_',yr,'_world.tif'), gdal=c('COMPRESS=NONE', 'TFW=YES'), datatype='INT1U')
  return('Done.\n')
})

p95_idx <- index |> purrr::map(1) |> terra::rast()
p95_whn <- index |> purrr::map(2) |> terra::rast()

# Two-years
yrs <- 1981:2023
p95_index <- lapply(X = yrs, FUN = function(yr){
  # Precipitation files filtered
  prc_flt <- prc[grep(pattern = c(yr,yr+1), x = prc)] # To implement
  rnf <- terra::rast(prc_flt)
  rsl <- terra::rapp(x = rnf, first = s1_ini_5km, last = s1_end_5km, fun = p95)
  terra::writeRaster(x = rsl, filename = paste0(root,'/yr_',yr,'_world.tif'), gdal=c('COMPRESS=NONE', 'TFW=YES'), datatype='INT1U')
  return('Done.\n')
})

# ------------------------------------------ #
# Two seasons
# ------------------------------------------ #
