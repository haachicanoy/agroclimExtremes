## ------------------------------------------ ##
## Resampling agro-climatic indices at 25 km
## By: Harold Achicanoy
## WUR & ABC
## Dec. 2023
## ------------------------------------------ ##

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,lubridate))

root <- '//CATALOGUE/WFP_ClimateRiskPr1'

# Raster template at 25 km
tmp <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5_25km.tif')

# Resampling function
resampling_indices <- function(index = 'spei-6'){
  
  cat(paste0('>>> Resampling index: ',index,'... \n'))
  cat(paste0('> One growing season\n'))
  # List files
  fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_10km'), pattern = paste0('^one_s1_',index,'_[0-9]*.*.tif$'), full.names = T)
  # Load indices at 10 km
  idx <- terra::rast(fls)
  # Resample indices at 25 km
  idx_25km <- terra::resample(x = idx, y = tmp, method = 'max')
  outfile <- paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_25km/one_s1_',index,'_25km.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = idx_25km, filename = outfile, overwrite = T)
  
  cat(paste0('> Two growing seasons: season 1.\n'))
  # List files
  fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_10km'), pattern = paste0('^two_s1_',index,'_[0-9]*.*.tif$'), full.names = T)
  # Load indices at 10 km
  idx <- terra::rast(fls)
  # Resample indices at 25 km
  idx_25km <- terra::resample(x = idx, y = tmp, method = 'max')
  outfile <- paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_25km/two_s1_',index,'_25km.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = idx_25km, filename = outfile, overwrite = T)
  
  cat(paste0('> Two growing seasons: season 2.\n'))
  # List files
  fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_10km'), pattern = paste0('^two_s2_',index,'_[0-9]*.*.tif$'), full.names = T)
  # Load indices at 10 km
  idx <- terra::rast(fls)
  # Resample indices at 25 km
  idx_25km <- terra::resample(x = idx, y = tmp, method = 'max')
  outfile <- paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_25km/two_s2_',index,'_25km.tif')
  dir.create(dirname(outfile),F,T)
  terra::writeRaster(x = idx_25km, filename = outfile, overwrite = T)
  cat('Done!\n\n')
  
}

# ------------------------------------------ #

resampling_indices(index = 'cdd')    # Maximum number of consecutive dry days
resampling_indices(index = 'spei-6') # SPEI-6

# ------------------------------------------ #














## GRAPHS TO CHECK

# Load clustering results
# fmc <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/clusters/WORLD_cdd_fmado_trimmed_clusters_100.tif'))
fmc <- terra::rast('D:/WORLD_fmado_trimmed_100_gs1.tif')
cdd_25km <- terra::rast(paste0(root,'/agroclimExtremes/agex_indices/agex_cdd/agex_cdd_25km/one_s1_cdd_25km.tif'))

# Combined data
dta <- c(cdd_25km, fmc)

# Combined data as data frame
dfm <- terra::as.data.frame(x = dta, xy = T, cell = T, na.rm = T)
names(dfm)[ncol(dfm)] <- 'cluster'

s_dfm_001 <- dfm[dfm$cluster == 1,]

s_dfm_17 <- dfm[dfm$cluster == 17,]

plot(x = 1979:2022, y = as.numeric(s_dfm_001[1,4:(ncol(s_dfm_001)-1)]), ty = 'l', ylim = c(0,300))
pacman::p_load(scales)
for(i in 2:nrow(s_dfm_001)){
  lines(x = 1979:2022, y = as.numeric(s_dfm_001[i,4:(ncol(s_dfm_001)-1)]), col = alpha('red', 0.05))
}; rm(i)

s_dfm_87 <- dfm[dfm$cluster == 87,]
for(i in 1:nrow(s_dfm_87)){
  lines(x = 1979:2022, y = s_dfm_87[i,4:(ncol(s_dfm_87)-1)], col = alpha('blue', 0.05))
}; rm(i)

s_dfm_36 <- dfm[dfm$cluster == 36,]
for(i in 1:nrow(s_dfm_36)){
  lines(x = 1979:2022, y = s_dfm_36[i,4:(ncol(s_dfm_36)-1)], col = alpha('forestgreen', 0.05))
}; rm(i)

r_17 <- fmc
r_17[r_17 != 17] <- NA
plot(terra::focal(x = r_17, w = 3, fun = mean, na.rm = F))

cdd <- terra::rast("//CATALOGUE/WFP_ClimateRiskPr1/agroclimExtremes/agex_indices/agex_cdd/agex_cdd_10km/one_s1_cdd_2022.tif")
## Graph of one index-year
wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
wrl <- terra::aggregate(wrl)
# my.palette <- RColorBrewer::brewer.pal(n = 20, name = 'Set1')
my.palette <- MetBrewer::met.brewer(name = 'Tam', n = 20)
tmp <- cdd[[terra::nlyr(cdd)]]
tmp <- terra::trim(x = tmp)
hist(terra::values(tmp))
qnt <- quantile(x = terra::values(tmp), probs = 0.98, na.rm = T)
terra::values(tmp)[terra::values(tmp) > qnt] <- qnt
png(filename = "D:/test.png", width = 3132, height = 2359, units = 'px')
plot(wrl, ext = terra::ext(tmp))
plot(tmp, add = T, col = my.palette, plg = list(cex = 5), cex.main = 7)
dev.off()

## Graph of global cluster
wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
wrl <- terra::aggregate(wrl)
# my.palette <- RColorBrewer::brewer.pal(n = 20, name = 'Set1')
my.palette <- MetBrewer::met.brewer(name = 'Signac', n = 100)
set.seed(1235); my.palette <- sample(x = my.palette, size = 100, replace = F)
tmp <- fmc
tmp <- terra::trim(x = tmp)
hist(terra::values(tmp))
png(filename = "D:/global_cluster.png", width = 3132, height = 2359, units = 'px')
plot(wrl, ext = terra::ext(tmp))
plot(tmp, add = T, col = my.palette, plg = list(cex = 5), cex.main = 7)
dev.off()


plot(ngs == 1)

gs1 <- ngs
gs1[gs1 != 1] <- NA
gs2 <- ngs
gs2[gs2 != 2] <- NA
gs2[!is.na(gs2)] <- 1

plot(terra::mask(x = fmc, mask = gs1))
terra::writeRaster(x = terra::mask(x = fmc, mask = gs1), filename = 'D:/WORLD_fmado_trimmed_100_gs1.tif', overwrite = T)

r <- terra::rast('D:/WORLD_fmado_trimmed_100_gs1.tif')
