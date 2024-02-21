# value of production vs agricultural gdp
library(terra)
library(scales)

vop <- terra::rast('C:/Users/haachicanoy/Downloads/spam2010V2r0_global_V_agg_VP_CROP_A.tif')

tmp <- vop
tmp[!is.na(tmp)] <- 1

gdp <- terra::rast('C:/Users/haachicanoy/Downloads/aggdp2010.tif')
gdp <- terra::resample(x = gdp, y = tmp, method = 'bilinear')
gdp <- terra::mask(x = gdp, mask = tmp)

dta <- c(vop, gdp)
dfm <- terra::as.data.frame(x = dta, xy = T, cell = T, na.rm = T)
plot(log(dfm[,4:5]), pch = 20, col = scales::alpha('black', 0.1))

tmp <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5_25km.tif')

vop_25km <- terra::resample(x = vop, y = tmp, method = 'bilinear')

# Spatial dependent clusters
fmd <- terra::rast('D:/WORLD_fmado_trimmed_100_gs1.tif')
plot(terra::zonal(x = vop_25km, z = fmd, fun = 'mean', as.raster = T, na.rm = T))

