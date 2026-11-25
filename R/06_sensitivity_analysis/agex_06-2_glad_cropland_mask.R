# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Sensitivity analysis: physical cropland masks (GLAD) before
# clustering (Reviewer 1)
# By: Harold Achicanoy
# WUR & ABC
# Created in November 2026
# --------------------------------------------------------------- #
# Steps:
# 1. Aggregate GLAD cropland fraction (Potapov et al. 2022, Nature Food)
#    from the 0.025 deg overview layers to the 0.25 deg grid as a MEAN
#    FRACTION (not nearest neighbour, which would keep the problem binary).
# 2. For each threshold (1, 5, 10, 20%) drop pixels below it and re-run
#    the F-madogram clustering for one and two growing seasons.
#    A reference run with no threshold (gladref) goes through the same
#    simplified pipeline (no outlier cleaning), so comparisons are
#    like-for-like.
# 3. Report how the number, size and geometry of clusters change, with
#    an explicit Central Africa (Congo Basin) summary.
# 4. Identify hotspots on each masked partition (crop diversity is
#    recomputed only over the retained pixels).
# 5. (R3.4) Cropland change 2003-2019 per cluster to flag hotspots located
#    in areas of rapid agricultural transformation.
# --------------------------------------------------------------- #

# R options and user-defined functions ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra, dplyr, cluster, future, furrr, landscapemetrics, aricode))
source('https://raw.githubusercontent.com/haachicanoy/agroclimExtremes/main/R/06_sensitivity/agex_06-0_sensitivity_functions.R')

# Optimal number of clusters (same criterion as 03-2: PAM + average silhouette
# on a 3% sample of pixels). k_range is an argument so it can be adapted to the
# smaller number of pixels retained after masking.
select_k <- function(x, k_range = 200:300, prc = 0.03, workers = 40, seed = 1245){
  set.seed(seed)
  smp <- sample(x = 1:ncol(x), size = max(2*max(k_range), round(prc * ncol(x))), replace = F)
  dd  <- get_fmado_dist(x[,smp]) |> cap_fmado_dist(); gc(T)
  future::plan(cluster, workers = workers, gc = T)
  sil_width <- furrr::future_map(.x = k_range, .f = function(k){
    model <- cluster::pam(x = dd, k = k)
    data.frame(k = k, avg_silhouette = model$silinfo$avg.width)
  }) |> dplyr::bind_rows()
  future:::ClusterRegistry('stop'); gc(verbose = F, reset = T, full = T)
  return(sil_width$k[which.max(sil_width$avg_silhouette)])
}

# Key arguments ----
root      <- '//CATALOGUE/AgroclimExtremes'
out_dir   <- paste0(root,'/agex_results/agex_sensitivity')
cls_dir   <- paste0(out_dir,'/clusters')
index     <- 'spei-6'
agr_dts   <- 'spam2010'
thrs      <- c(0, 0.01, 0.05, 0.10, 0.20) # 0 = reference without mask
k_mode    <- 'optimal'                     # 'optimal' (silhouette) or 'baseline' (same K as main run)
glad_dir  <- paste0(root,'/agex_raw_data/glad_cropland') # CHECK: GLAD 0.025 deg overview layers (percent cropland)
glad_yrs  <- c(2003, 2019)
cafr_ext  <- terra::ext(c(12, 31, -12, 6))   # Congo Basin / Central Africa box (lon/lat)
dir.create(cls_dir, F, T)

# Template 25 km grid
tmp_25km <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_',index,'_combined_fmadogram_clean.tif'))

# 1. GLAD cropland fraction at 0.25 deg ----
for(yr in glad_yrs){
  outfile <- paste0(root,'/agex_results/agex_vulnerability/glad_cropland_fraction_',yr,'_25km.tif')
  if(!file.exists(outfile)){
    glad <- terra::rast(list.files(path = glad_dir, pattern = paste0(yr,'.*\\.tif$'), full.names = T)[1]) # CHECK file pattern
    glad <- terra::classify(glad, rcl = cbind(101, Inf, NA)) # remove no-data codes, if any
    glad_25km <- terra::aggregate(x = glad, fact = round(terra::res(tmp_25km)[1]/terra::res(glad)[1]), fun = 'mean', na.rm = T)
    glad_25km <- terra::resample(x = glad_25km, y = tmp_25km, method = 'near') / 100
    names(glad_25km) <- 'cropland_fraction'
    terra::writeRaster(x = glad_25km, filename = outfile, overwrite = T)
  }
}; rm(yr, outfile)
glad_frc <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/glad_cropland_fraction_',max(glad_yrs),'_25km.tif'))

# 2. Re-run clustering with cropland thresholds ----
stp <- data.frame(gs = c('one','two','two'), season = c(1,1,2))

for(thr in thrs){
  scn_lbl <- ifelse(thr == 0, 'gladref', sprintf('glad%02d', round(100*thr)))
  outfile <- paste0(cls_dir,'/agex_global_',index,'_combined_fmadogram_',scn_lbl,'.tif')
  if(file.exists(outfile)) next
  cat(paste0('Processing threshold: ',thr,'\n'))

  cls_set <- lapply(1:nrow(stp), function(i){
    idx_roi <- load_index_ts(root, index, stp$gs[i], stp$season[i])
    idx_roi_fnl <- terra::app(x = idx_roi, fun = function(i, ff) ff(i), cores = 20, ff = transformations)
    idx_roi_ntrd <- terra::as.data.frame(x = idx_roi_fnl, xy = T, cell = T, na.rm = T); rm(idx_roi_fnl); gc(T)
    # Cropland mask
    frc <- terra::extract(x = glad_frc, y = idx_roi_ntrd$cell)[,1]
    keep <- !is.na(frc) & frc >= thr
    if(thr == 0) keep <- rep(TRUE, nrow(idx_roi_ntrd))
    idx_roi_ntrd <- idx_roi_ntrd[keep,]
    x <- t(idx_roi_ntrd[,4:ncol(idx_roi_ntrd)]); gc(T)
    k <- switch(k_mode,
                optimal  = select_k(x),
                baseline = get_baseline_k(root, index, stp$gs[i], stp$season[i]))
    cls <- fmado_cluster(x, k)
    r <- idx_roi[[1]]; terra::values(r) <- NA; r[idx_roi_ntrd$cell] <- cls
    names(r) <- 'extreme_cluster'
    terra::writeRaster(x = r, filename = paste0(cls_dir,'/agex_global_',index,'_',stp$gs[i],'_s',stp$season[i],'_fmadogram_k',k,'_',scn_lbl,'.tif'), overwrite = T)
    return(r)
  })

  # Combine two growing seasons (as 03-4, without outlier cleaning):
  # keep pixels clustered in both seasons; drop combinations smaller than the
  # smallest one-season cluster.
  lwr_lmt <- min(table(as.numeric(terra::values(cls_set[[1]], na.rm = T))))
  dfm <- terra::as.data.frame(x = c(cls_set[[2]], cls_set[[3]]), cell = T, na.rm = T)
  names(dfm)[2:3] <- c('s1','s2')
  dfm$combined <- paste0(dfm$s1,'-',dfm$s2)
  aux <- base::as.data.frame(sort(table(dfm$combined), decreasing = T)); names(aux) <- c('combined','count')
  aux <- aux[aux$count >= lwr_lmt,]
  aux$extreme_cluster <- 1:nrow(aux)
  dfm <- dplyr::inner_join(x = dfm, y = aux[,c('combined','extreme_cluster')], by = 'combined')

  # Merge one and two growing seasons (as 03-5) with non-overlapping ids
  agex_sgn <- cls_set[[1]]
  n_gs1 <- max(terra::values(agex_sgn, na.rm = T))
  agex_sgn[dfm$cell] <- dfm$extreme_cluster + n_gs1
  # Re-number consecutively
  vls <- sort(unique(terra::values(agex_sgn, na.rm = T)))
  agex_sgn <- terra::classify(agex_sgn, rcl = cbind(vls, 1:length(vls)))
  names(agex_sgn) <- 'extreme_cluster'
  terra::writeRaster(x = agex_sgn, filename = outfile, overwrite = T)
  rm(cls_set, dfm, aux, agex_sgn, vls, n_gs1, lwr_lmt)
}; rm(thr, outfile, scn_lbl)

# 3. Cluster geometry per threshold ----
cls_fls <- list.files(path = cls_dir, pattern = paste0('agex_global_',index,'_combined_fmadogram_glad.*\\.tif$'), full.names = T)
ref <- terra::rast(grep('gladref', cls_fls, value = T))
geo <- lapply(cls_fls, function(f){
  r <- terra::rast(f); names(r) <- 'extreme_cluster'
  are <- terra::zonal(x = terra::cellSize(r, unit = 'km') |> terra::mask(r), z = r, fun = 'sum', na.rm = T)
  npt <- landscapemetrics::lsm_c_np(landscape = r) # patches per cluster: fragmentation
  r_ca <- terra::crop(r, cafr_ext)
  are_ca <- terra::zonal(x = terra::cellSize(r_ca, unit = 'km') |> terra::mask(r_ca), z = r_ca, fun = 'sum', na.rm = T)
  dfm <- terra::as.data.frame(c(ref, r), cell = T, na.rm = T)
  data.frame(scenario            = gsub('.*_(glad[a-z0-9]+)\\.tif$', '\\1', f),
             n_clusters          = nrow(are),
             n_pixels            = sum(!is.na(terra::values(r))),
             total_area_km2      = sum(are[,2]),
             median_area_km2     = stats::median(are[,2]),
             max_area_km2        = max(are[,2]),
             mean_patches        = mean(npt$value),
             ARI_vs_ref          = aricode::ARI(dfm[,2], dfm[,3]),
             cafr_n_clusters     = nrow(are_ca),
             cafr_area_km2       = sum(are_ca[,2]),
             cafr_median_area_km2 = stats::median(are_ca[,2]))
}) |> dplyr::bind_rows()
utils::write.csv(x = geo, file = paste0(out_dir,'/agex_glad_cluster_geometry_',index,'.csv'), row.names = F)
geo

# 4. Hotspots on masked partitions ----
for(f in cls_fls){
  scn_lbl <- gsub('.*_(glad[a-z0-9]+)\\.tif$', '\\1', f)
  run_hotspots(root = root, index = index, agr_dts = agr_dts, cls_file = f,
               scenario = paste0(index,'_',agr_dts,'_',scn_lbl), out_dir = out_dir)
}; rm(f, scn_lbl)

# 5. Cropland change 2003-2019 within baseline clusters (R3.4) ----
frc_0 <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/glad_cropland_fraction_',min(glad_yrs),'_25km.tif'))
frc_1 <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/glad_cropland_fraction_',max(glad_yrs),'_25km.tif'))
dlt <- frc_1 - frc_0; names(dlt) <- 'cropland_fraction_change'
bsl <- utils::read.csv(paste0(out_dir,'/agex_hotspots_',index,'_',agr_dts,'_glad00.csv'))
cls <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_',index,'_combined_fmadogram_clean.tif')); names(cls) <- 'extreme_cluster'
chg <- terra::zonal(x = c(dlt, abs(dlt) |> stats::setNames('abs_change')), z = cls, fun = 'mean', na.rm = T)
chg$rapid_change <- chg$abs_change >= stats::quantile(chg$abs_change, .90, na.rm = T) # top decile of clusters
bsl <- dplyr::left_join(bsl, chg, by = 'extreme_cluster')
utils::write.csv(x = bsl, file = paste0(out_dir,'/agex_hotspots_',index,'_',agr_dts,'_cropland_change.csv'), row.names = F)
# Hotspots located in rapidly changing agricultural areas
bsl |>
  dplyr::filter(food_production_hotspot | economic_value_hotspot) |>
  dplyr::summarise(n_hotspots = dplyr::n(), n_rapid_change = sum(rapid_change, na.rm = T))
