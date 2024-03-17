## ------------------------------------------ ##
## Thermal Humidity Index (multiple species)
## By: Harold Achicanoy
## WUR & ABC
## Mar. 2024
## ------------------------------------------ ##

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

# Source agro-climatic indices functions
source('https://raw.githubusercontent.com/haachicanoy/agroclimExtremes/main/agex_00_indices.R')
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')

# Root directory
root <- '//CATALOGUE/WFP_ClimateRiskPr1'

# Calculate THI function
calc_thi <- function(yr, mn){
  # Output files
  outfile1 <- paste0(out_dir,'/THI_mean-',yr,'-',mn,'.tif')
  outfile2 <- paste0(out_dir,'/THI_max-',yr,'-',mn,'.tif')
  outfile3 <- paste0(out_dir,'/daily/THI_daily-',yr,'-',mn,'.tif')
  file.exists(c(outfile1,outfile2,outfile3))
  cat(outfile2,'\n')
  if(sum(file.exists(c(outfile1,outfile2,outfile3))) < 3){
    dir.create(dirname(outfile3),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    
    # Files
    tav_fls <- paste0(ae5tx_pth,'/Temperature-Air-2m-Max-24h_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.1.nc')
    tav_fls <- tav_fls[file.exists(tav_fls)]
    rhm_fls <- paste0(ae5rh_pth,'/Relative-Humidity-2m-12h_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.1.nc')
    rhm_fls <- rhm_fls[file.exists(rhm_fls)]
    
    # Read variables
    tav <- terra::rast(tav_fls)
    tav <- tav - 273.15
    rhm <- terra::rast(rhm_fls)
    
    # Calculate thermal humidity index for multiple species
    THI <- terra::lapp(x = terra::sds(tav, rhm), fun = thi)
    terra::time(THI) <- dts
    THI_avg <- mean(THI, na.rm = T)
    THI_max <- max(THI, na.rm = T)
    
    # Write output
    terra::writeRaster(THI_avg, outfile1)
    terra::writeRaster(THI_max, outfile2)
    terra::writeRaster(THI, outfile3)
    
    # Clean up
    rm(tmx, rhm, THI, THI_avg, THI_max)
    gc(verbose=F, full=T, reset=T)
  }
}

#input data paths for agera5
ae5tx_pth <- paste0(root,'/1.Data/AgERA5/2m_temperature-24_hour_mean') # Mean temperature
ae5rh_pth <- paste0(root,'/1.Data/AgERA5/2m_relative_humidity') # Relative humidity

#output path
out_dir <- paste0(root,'/agroclimExtremes/agex_raw_data/monthly_thi')

# Historical setup
yrs <- 1979:2023
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp |>
  dplyr::arrange(yrs, mns) |>
  base::as.data.frame()

1:nrow(stp) |>
  purrr::map(.f = function(i){
    calc_thi(yr = stp$yrs[i], mn = stp$mns[i])
    gc(verbose=F, full=T, reset=T)
    terra::tmpFiles(remove = T)
  })
