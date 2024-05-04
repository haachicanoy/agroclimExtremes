# ------------------------------------------ #
# Monthly water balance at 10 km
# By: Harold Achicanoy
# WUR & ABC
# Feb. 2024
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(tidyverse,terra))

# Root directory
root <- '//CATALOGUE/WFP_ClimateRiskPr1'

# Historical setup
yrs <- 1979:2023
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>%
  dplyr::arrange(yrs, mns) %>%
  base::as.data.frame()
# Input directories
evp_pth <- paste0(root,'/agroclimExtremes/agex_raw_data/monthly_evapotranspiration')
prc_pth <- paste0(root,'/agroclimExtremes/agex_raw_data/monthly_precipitation')
# Output directory
out_pth <- paste0(root,'/agroclimExtremes/agex_raw_data/monthly_balance')

calc_balance <- function(yr, mn){
  evp <- terra::rast(paste0(evp_pth,'/ET-',yr,'-',mn,'.tif'))
  prc <- terra::rast(paste0(prc_pth,'/prec-',yr,'-',mn,'.tif'))
  bal <- prc - evp
  terra::writeRaster(x = bal, filename = paste0(out_pth,'/bal-',yr,'-',mn,'.tif'))
}

# loop for each year and month
1:nrow(stp) |>
  purrr::map(.f = function(i){calc_balance(yr = stp$yrs[i], mn = stp$mns[i]); gc(T)})
