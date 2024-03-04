## Priestley-Taylor Model Evapotranspiration (ETP)
## By: Cesar Saavedra
## Junio, 2023

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))
#
root <- '//CATALOGUE/WFP_ClimateRiskPr1'
cnty <- terra::vect(paste0(root,'/1.Data/shps/GTM/GTM_GADM1.shp'))

# Calculate FAO-Penman-Montieth ET function
calc_ET0 <- function(yr, mn){
  outfile1 <- paste0(out_dir,'/Priestley-Taylor/ET_0_monthly/ET-',yr,'-',mn,'.tif')
  # outfile2 <- paste0(out_dir,'/Priestley-Taylor/ET_0_daily/',yr,'/ET-',yr,'-',mn,'-')
  cat(outfile1, "\n")
   if(!file.exists(outfile1)){
    dir.create(dirname(outfile1),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    if(as.numeric(yr) > 2020 & mn == '02'){
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-28')), by = 'day')
    } else {
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    }
    # Files
    tx_fls <- paste0(tx_pth,'Temperature-Air-2m-Max-24h_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0','.nc')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    tm_fls <- paste0(tm_pth,'Temperature-Air-2m-Min-24h_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0','.nc')
    tm_fls <- tm_fls[file.exists(tm_fls)]
    sr_fls <- paste0(sr_pth,'Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0','.nc')
    sr_fls <- sr_fls[file.exists(sr_fls)]
    #
    ref <- terra::rast(paste0(root,"/1.Data/chirps-v2.0.1981.01.1.tif"))
    # Read variables
    #
    tmx <- terra::rast(tx_fls)
    tmx <- terra::resample(tmx,ref) %>% terra::crop(terra::ext(cnty)) %>% terra::mask(cnty)
    tmx[tmx == -9999] <- NA
    tmx <- tmx - 273.15
    #
    tmn <- terra::rast(tm_fls)
    tmn <- terra::resample(tmn,ref) %>% terra::crop(terra::ext(cnty)) %>% terra::mask(cnty)
    tmn[tmn == -9999] <- NA
    tmn <- tmn - 273.15
    #
    sr <- terra::rast(sr_fls)
    sr <- terra::resample(sr,ref) %>% terra::crop(terra::ext(cnty)) %>% terra::mask(cnty)
    sr <- sr/1000000 # solar radiation from j to kj
    sr <- (1/2.45)*(sr) # Solar radiation from kj to mm day-1
    # 
    tmean <- (tmx+tmn)/2
    #
    elev <- terra::rast(paste0(root,'/1.Data/DEM_10km_Guatemala.tif'))
    lst <- list()
    for(i in 1:length(dts)){
      print(i)
      lst[[i]] <- elev
    }
    elv <- terra::rast(lst)
    elv <- resample(elv,ref) %>% terra::crop(terra::ext(cnty)) %>% terra::mask(cnty)
    # calculate ET constants 
    # a = 1.26                      # Priestley and Taylor evaporative coefficient
    a <- tmean
    terra::values(a) <- 1.26  
    #
    # b = 2.45                      # Latent heat of vaporization at 20°C 
    b <- tmean
    terra::values(b) <- 2.45
    #
    d = (4098*(0.6108*(exp((17.27*tmean)/(tmean+237.3)))))/((tmean+237.3)*(tmean+237.3)) # Relationship between saturation vapor pressure and air temperature
    #
    e = (-0.000007*elv)+0.0666    # Psychrometric constant
    # Evapotranspiration function
    ET_0_idx <- function(a, b, d, e, sr){
      ET_pt = a*(d/(d+e))*(sr/b)  #Priestley-Taylor Model
      return(ET_pt)
    }
    # Calculate evapotranspiration 
    ET_0 <- terra::lapp(x = terra::sds(a, b, d, e, sr), 
                        fun = ET_0_idx) %>% sum()
    # Files names
    # fn <- paste0(1:nlyr(ET_0), '.tif') 
    terra::writeRaster(ET_0, outfile1, overwrite=T)
    # terra::writeRaster(ET_0, paste0(outfile2,fn), overwrite=T)
   }
  cat("Done", "\n", ":)", "\n")
}

# Historical setup
yrs <- 1981:2022
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>%
  dplyr::arrange(yrs, mns) %>%
  base::as.data.frame()
tx_pth <- paste0(root,'/1.Data/ERA5/2m_temperature-24_hour_maximum/') # Maximum temperature
tm_pth <- paste0(root,'/1.Data/ERA5/2m_temperature-24_hour_minimum/') # Minimun temperature
sr_pth <- paste0(root,'/1.Data/ERA5/solar_radiation_flux/')           # Solar radiation

out_dir <- '//CATALOGUE/AgriLACRes_WP2/1.Data/Guatemala/SPEI_Historical_analysis'
1:nrow(stp) %>%
  purrr::map(.f = function(i){calc_ET0(yr = stp$yrs[i], mn = stp$mns[i]); gc(verbose=F, full=T, reset=T)})
