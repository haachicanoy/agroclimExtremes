library(terra)

root <- '//CATALOGUE/AgroclimExtremes'
cdd <- terra::rast(paste0(root,'/agex_indices/agex_cdd/agex_cdd_25km/one_s1_cdd_25km.tif'))
cdd <- cdd[[2:(terra::nlyr(cdd))]]
cdd[is.infinite(cdd)] <- 0

spei <- terra::rast(paste0(root,'/agex_indices/agex_spei-6/agex_spei-6_25km/one_s1_spei-6_25km.tif'))

idx_cor <- terra::lapp(x = terra::sds(cdd,spei), fun = cor, method = 'spearman')
plot(idx_cor)

terra::writeRaster(x = idx_cor, filename = 'D:/correlation__cdd_vs_spei.tif')

## Load Palmer Drought Index
palmer <- terra::rast('C:/Users/haachicanoy/Downloads/palmer_index/scPDSI.cru_ts4.07early1.1901.2022.cal_1901_22.bams.2023.GLOBAL.IGBP.WHC.1901.2022.nc')
palmer <- palmer[[terra::time(palmer) > as.Date('1979-12-01')]]

palmer <- terra::tapp(x = palmer, index = lubridate::year(terra::time(palmer)), fun = mean, cores = 20)
palmer <- terra::resample(x = palmer, y = spei, method = 'bilinear')
palmer <- terra::mask(x = palmer, mask = spei[[1]])

idx_cor2 <- terra::lapp(x = terra::sds(palmer,spei), fun = cor, method = 'spearman')
plot(idx_cor2)

terra::writeRaster(x = idx_cor2, filename = 'D:/correlation__palmeravg_vs_spei.tif')

## Load THI max
thi <- terra::rast(list.files(path = paste0(root,'/agex_raw_data/monthly_thi'), pattern = 'THI_max', full.names = T))
nms <- basename(list.files(path = paste0(root,'/agex_raw_data/monthly_thi'), pattern = 'THI_max', full.names = T))
nms <- gsub('THI_max-','',nms)
nms <- gsub('.tif','',nms)
dts <- paste0(nms,'-01')
dts <- as.Date(dts)

thi <- terra::tapp(x = thi, index = lubridate::year(dts), fun = max, cores = 20)
thi <- terra::resample(x = thi, y = spei[[1]], method = 'bilinear', threads = T)
thi <- terra::mask(x = thi, mask = spei[[1]])
thi <- thi[[2:(terra::nlyr(thi)-1)]]

idx_cor3 <- terra::lapp(x = terra::sds(thi,spei), fun = cor, method = 'spearman')
plot(idx_cor3)

terra::writeRaster(x = idx_cor3, filename = 'D:/correlation__thi_vs_spei.tif')


# ---------------------------------------------------------- # 
# Coordinate 1
# ---------------------------------------------------------- # 

cll1 <- terra::cellFromRowCol(cdd , row = 512, col = 1305)
crd1 <- terra::xyFromCell(spei, cll1)

thi_vct1 <- terra::extract(thi, crd1) |> as.numeric()
spei_vct1 <- terra::extract(spei, crd1) |> as.numeric()

cor(x = thi_vct1, y = spei_vct1, method = 'spearman')

plot(scale(thi_vct1), ty = 'l', ylim = c(-3,3))
lines(spei_vct1, col = 'red')
abline(h = 1.5, col = 'blue')

plot(thi_vct1, spei_vct1, pch = 20)
abline(h = 1.5, col = 'blue')

# Bivariate empirical distribution
bed1 <- robusTest::ecdf2D(x = thi_vct1, y = spei_vct1)
plot(terra::rast(bed1$ecdf))

# ---------------------------------------------------------- # 
# Coordinate 2
# ---------------------------------------------------------- # 

cll2 <- terra::cellFromRowCol(cdd , row = 360, col = 1190)
crd2 <- terra::xyFromCell(cdd, cll2)

cdd_vct2 <- terra::extract(cdd, crd2) |> as.numeric()
spei_vct2 <- terra::extract(spei, crd2) |> as.numeric()

cor(x = cdd_vct2, y = spei_vct2, method = 'spearman')

plot(scale(cdd_vct2), ty = 'l')
lines(spei_vct2, col = 'red')
abline(h = 1.5, col = 'blue')

plot(cdd_vct2, spei_vct2, pch = 20)
abline(h = 1.5, col = 'blue')

# Bivariate empirical distribution
bed2 <- robusTest::ecdf2D(x = cdd_vct2, y = spei_vct2)
plot(terra::rast(bed2$ecdf))

plot(terra::rast(bed1$ecdf))
plot(terra::rast(bed2$ecdf))
