# ------------------------------------------ #
# FAO-Penman-Monteith Evapotranspiration
# By: Cesar Saavedra & Harold Achicanoy
# ABC
# Dec. 2023
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
rm(list = ls()); gc(TRUE)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- '//CATALOGUE/WFP_ClimateRiskPr1'
cnty <- terra::vect(paste0(root,'/1.Data/shps/GTM/GTM_GADM1.shp'))

# Calculate FAO-Penman-Monteith ET function
calc_ET0 <- function(yr, mn){
  outfile <- paste0(out_dir,'/FAO-Penman-Monteith/ET_0_monthly/ET-',yr,'-',mn,'.tif')
  # outfile2 <- paste0(out_dir,'/FAO-Penman-Monteith/ET_0_daily/',yr,'/ET-',yr,'-',mn,'-')
  cat(outfile, "\n")
  if(!file.exists(outfile)){
    dir.create(dirname(outfile),F,T)
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
    tav_fls <- paste0(tav_pth,'Temperature-Air-2m-Mean-24h_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0','.nc')
    tav_fls <- tav_fls[file.exists(tav_fls)]
    dp_fls <- paste0(dp_pth,'Dew-Point-Temperature-2m-Mean_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0','.nc')
    dp_fls <- dp_fls[file.exists(dp_fls)]
    sr_fls <- paste0(sr_pth,'Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0','.nc')
    sr_fls <- sr_fls[file.exists(sr_fls)]
    ws_fls <- paste0(ws_pth,'Wind-Speed-10m-Mean_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0','.nc')
    ws_fls <- ws_fls[file.exists(ws_fls)]
    #
    ref <- terra::rast(paste0(root,"/1.Data/chirps-v2.0.1981.01.1.tif"))
    # Read variables
    #
    tmx <- terra::rast(tx_fls)
    tmx <- terra::resample(tmx,ref) %>%  terra::crop(terra::ext(cnty)) %>% terra::mask(cnty)
    tmx[tmx == -9999] <- NA
    tmx <- tmx - 273.15
    #
    tmn <- terra::rast(tm_fls)
    tmn <- terra::resample(tmn,ref) %>% terra::crop(terra::ext(cnty)) %>% terra::mask(cnty)
    tmn[tmn == -9999] <- NA
    tmn <- tmn - 273.15
    #
    tav <- terra::rast(tav_fls)
    tav <- terra::resample(tav,ref) %>% terra::crop(terra::ext(cnty)) %>% terra::mask(cnty)
    tav[tav == -9999] <- NA
    tav <- tav - 273.15
    #
    dp <- terra::rast(dp_fls)
    dp <- terra::resample(dp,ref) %>% terra::crop(terra::ext(cnty)) %>% terra::mask(cnty)
    dp[dp == -9999] <- NA
    dp <- dp - 273.15
    #
    sr <- terra::rast(sr_fls)
    sr <- terra::resample(sr,ref) %>% terra::crop(terra::ext(cnty)) %>% terra::mask(cnty)
    sr <- sr/1000000 # solar radiation from j to kj
    sr <- 1/2.45*(sr) # Solar radiation from kj to mm day-1
    #
    ws <- terra::rast(ws_fls)
    ws <- terra::resample(ws,ref) %>% terra::crop(terra::ext(cnty)) %>% terra::mask(cnty)
    ws_2m <- ws*((4.87)/(log((67.8*10)-5.42)))
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
    p <- (101.3)*((293-0.0065*elv)/(293))^5.26   # Presion
    cps <- ((1.013*10^-3)*p)/(0.622*2.45)        # Constante psicrometrica
    etmx <- 0.6108*exp((17.27*tmx)/(tmx+237.3))  # Presion media de vapor de la saturacion
    etmn <- 0.6108*exp((17.27*tmn)/(tmn+237.3))  # Presion media de vapor de la saturacion
    es <- (etmx+etmn)/2                          # Presion media de vapor de la saturacion
    Dlt <- (4098*(0.6108*(exp((17.27*tmean)/(tmean+237.3)))))/((tmean+237.3)^2) # Delta/Pendiente de la curva de presion de saturacion de vapor
    ea <- 0.6108*(exp((17.27*dp)/(dp+237.3)))    # Presion real de vapor derivada del punto de rocio
    # Evapotranspiration function
    ET_0_idx <- function(Dlt, sr, cps, tmean, ws_2m, es){
      et_0 <- ((0.408*Dlt*sr)+(cps*(900/(tmean+273))*ws_2m*es))/(Dlt+(cps*(1+(0.34*ws_2m))))
      return(et_0)
    }
    # Calculate evapotranspiration 
    ET_0 <- terra::lapp(x = terra::sds(Dlt, sr, cps, tmean, ws_2m, es), 
                        fun = ET_0_idx) %>% sum()
    # Files names
    # fn <- paste0(1:nlyr(ET_0), '.tif')
    terra::writeRaster(ET_0, outfile, overwrite=T)
    # terra::writeRaster(ET_0, paste0(outfile2,fn), overwrite=T)
  }
  cat("Done", "\n", ":)", "\n")
}

# Historical setup
yrs <- 1981:2022
mns <- c(paste0('0',1:9),10:12)
# dys <- c(paste0('0'1:9),10:30)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>%
  dplyr::arrange(yrs, mns) %>%
  base::as.data.frame()
tx_pth <- paste0(root,'/1.Data/ERA5/2m_temperature-24_hour_maximum/') # Maximum temperature
tm_pth <- paste0(root,'/1.Data/ERA5/2m_temperature-24_hour_minimum/') # Minimun temperature
tav_pth <- paste0(root,'/1.Data/ERA5/2m_temperature-24_hour_mean/')   # Mean temperature
sr_pth <- paste0(root,'/1.Data/ERA5/solar_radiation_flux/')           # Solar radiation
ws_pth <- paste0(root,'/1.Data/ERA5/10m_wind_speed/')                 # Wind speed
dp_pth <- paste0(root,'/1.Data/ERA5/2m_dewpoint_temperature/')        # Dewpoint temperature

out_dir <- '//CATALOGUE/AgriLACRes_WP2/1.Data/Guatemala/SPEI_Historical_analysis'
# loop for each year and month
1:nrow(stp) %>%
  purrr::map(.f = function(i){calc_ET0(yr = stp$yrs[i], mn = stp$mns[i]); gc(verbose=F, full=T, reset=T)})
#
# SPEI::hargreaves()
