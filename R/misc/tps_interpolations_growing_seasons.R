set.seed(1235)
seeds <- round(runif(n = 100, min = 0, max = 100000))

aux <- n_seasons
terra::values(aux) <- NA

# Run the process over 100 random samples of 1000 coordinates
tsp_simulation <- lapply(seeds, function(seed){
  set.seed(seed)
  smp <- sample(x = n_seasons_dfm$cell, size = 1000, replace = F, prob = n_seasons_dfm$probability)
  sub_dfm <- n_seasons_dfm[n_seasons_dfm$cell %in% smp,]
  tps_fit <- fields::Tps(x = as.matrix(sub_dfm[,c('x','y')]), Y = sub_dfm$nseasons)
  aux <- terra::rast(n_seasons)
  aux <- terra::interpolate(aux, tps_fit)
  aux <- terra::mask(aux, n_seasons)
  return(aux)
})
tsp_simulation <- terra::rast(tsp_simulation)

tps_res <- terra::rast('//CATALOGUE/WFP_ClimateRiskPr1/agroclimExtremes/agex_raw_data/tps_simulation.tif')
tst <- terra::values(tps_res, na.rm = T)
colnames(tst) <- paste0('sim_',1:ncol(tst))
summary(tst)
min(tst)
max(tst)

hist(terra::values(tps_res[[1]], na.rm = T))

plot(ecdf(x = terra::values(tps_res[[3]], na.rm = T)))

terra::writeRaster(x = sum(tps_res <= 1.5), 'D:/one-gs_sites.tif')
terra::writeRaster(x = sum(tps_res > 1.5), 'D:/two-gs_sites.tif')