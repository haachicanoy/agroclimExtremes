# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Sensitivity analysis: shared functions
# By: Harold Achicanoy
# WUR & ABC
# Created in November 2026
# --------------------------------------------------------------- #
# This script only defines functions. It is sourced by:
# - agex_06-1_hotspots_identification.R
# - agex_06-2_glad_cropland_mask.R
# - agex_06-3_compare_scenarios.R
# - agex_06-4_nonstationarity.R
# Functions reproduce the logic of 03-1 (F-madogram), 04 (metrics)
# and 05 figure scripts (4x4 bivariate classification) but are
# parameterised by drought index, agricultural dataset and the
# extreme drought clusters raster, so every scenario follows the
# exact same protocol.
# --------------------------------------------------------------- #

suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra, dplyr, tidyr, purrr, Kendall, quantreg,
                                parallelDist, fastcluster, entropy, aricode,
                                rnaturalearth))

# --------------------------------------------------------------- #
# 1. Time series transformations and F-madogram (from 03-1)
# --------------------------------------------------------------- #

# Scale and, if the Mann-Kendall trend is significant, remove it with loess.
# detrend = FALSE keeps the scaled raw series (used in 06-4).
transformations <- function(x, detrend = TRUE){
  if(!all(is.na(x))){
    x_scl <- as.numeric(scale(x))
    if(detrend){
      trd <- Kendall::MannKendall(x_scl)$sl
      cnd <- ifelse(trd <= 0.05, T, F)
      if(cnd){
        dfm <- data.frame(x = 1:length(x_scl), y = x_scl)
        fit <- loess(y ~ x, data = dfm)
        x_ntrd <- dfm$y - fit$fitted
      } else {
        x_ntrd <- x_scl
      }
    } else {
      x_ntrd <- x_scl
    }
  } else {
    x_ntrd <- rep(NA, length(x))
  }
  return(x_ntrd)
}

# F-madogram distance (Bernard et al. 2013; Saunders et al. 2021)
get_fmado_dist <- function(x){
  Nnb <- ncol(x)
  Tnb <- nrow(x)
  V <- array(NaN, dim = c(Tnb,Nnb))
  for(p in 1:Nnb){
    x.vec <- as.vector(x[,p])
    if(all(is.na(x.vec))){
      V[,p] <- x.vec
      next
    }
    Femp <- ecdf(x.vec)(x.vec)
    V[,p] <- Femp
  }
  DD_fmado <- parallelDist::parDist(x = t(V), method = 'manhattan', diag = F, upper = F)/(2*Tnb)
  return(DD_fmado)
}

# Trim distribution of F-madogram distances greater than 1/6
cap_fmado_dist <- function(DD_fmado){
  DD_fmado_adj <- pmin(DD_fmado, rep(1/6, length(DD_fmado)))
  if(any(DD_fmado_adj > 1/6, na.rm = TRUE) == 1)
    stop("ERROR: cap_fmado_dist didn't work")
  return(DD_fmado_adj)
}

# Cluster a matrix of time series (rows = years, cols = pixels) with a fixed K.
# Returns a vector of cluster ids (one per pixel/column)
fmado_cluster <- function(x, k){
  dd  <- get_fmado_dist(x) |> cap_fmado_dist(); gc(T)
  hcl <- fastcluster::hclust(dd, method = 'ward.D2'); rm(dd); gc(T)
  cls <- stats::cutree(tree = hcl, k = k); rm(hcl); gc(T)
  return(cls)
}

# Load the yearly extreme SPEI layers of one growing-season setting
# (gs = 'one'/'two', season = 1/2) as produced by 02-3
load_index_ts <- function(root, index, gs, season, yrs = 1980:2022){
  fls <- list.files(path = paste0(root,'/agex_indices/agex_',index,'/agex_',index,'_25km'),
                    pattern = paste0(gs,'_s',season,'_',index,'_'), full.names = T)
  fls <- fls[grepl(pattern = paste0('_(',paste(yrs, collapse = '|'),').tif$'), x = fls)] |> sort()
  idx <- terra::rast(fls)
  return(idx)
}

# Number of clusters used in the baseline run of a growing-season setting
# (parsed from the file name written by 03-1: ..._fmadogram_k<K>.tif)
get_baseline_k <- function(root, index, gs, season){
  fl <- list.files(path = paste0(root,'/agex_results/agex_results_clusters'),
                   pattern = paste0('agex_global_',index,'_',gs,'_s',season,'_fmadogram_k[0-9]+.tif$'),
                   full.names = F)
  if(length(fl) == 0) stop(paste0('No baseline K found for ',index,' ',gs,' s',season))
  k <- as.numeric(gsub('.*_k([0-9]+).tif$', '\\1', fl[1]))
  return(k)
}

# --------------------------------------------------------------- #
# 2. Extreme drought intensification (from 04)
# --------------------------------------------------------------- #
# IMPORTANT: 02-3 stores -1 * min(SPEI) per growing season. Therefore the
# 95th percentile of the stored series is the 5th percentile of SPEI and a
# POSITIVE slope means intensifying droughts.
get_ext_trend <- function(x, tau = .95){
  if(!all(is.na(x))){
    x <- as.numeric(na.omit(x))
    dfm <- data.frame(time = 1:length(x), ts = x)
    qrf <- quantreg::rq(formula = ts ~ time, tau = tau, data = dfm)
    y <- as.numeric(qrf$coefficients[2])
  } else { y <- NA }
  return(y)
}

# Pixel-level extreme trend for a given index, merging one and two growing seasons.
# yrs allows computing the trend on sub-periods if needed.
calc_extreme_trend <- function(root, index, yrs = 1980:2022, cores = 20){
  sfx <- ifelse(identical(yrs, 1980:2022), '', paste0('_',min(yrs),'-',max(yrs)))
  outfile <- paste0(root,'/agex_results/agex_vulnerability/',toupper(index),'_extreme_trend',sfx,'.tif')
  if(!file.exists(outfile)){
    dir.create(path = dirname(outfile), F, T)
    stp <- data.frame(gs = c('one','two','two'), season = c(1,1,2))
    sxt <- lapply(1:nrow(stp), function(i){
      idx <- load_index_ts(root, index, stp$gs[i], stp$season[i], yrs)
      terra::app(x = idx, fun = function(i, ff) ff(i), cores = cores, ff = get_ext_trend)
    })
    idx_two_sxt <- mean(c(sxt[[2]], sxt[[3]]))
    idx_sxt <- terra::merge(sxt[[1]], idx_two_sxt); rm(sxt, idx_two_sxt)
    names(idx_sxt) <- 'extreme_trend'
    terra::writeRaster(x = idx_sxt, filename = outfile, overwrite = T)
  } else {
    idx_sxt <- terra::rast(outfile)
    names(idx_sxt) <- 'extreme_trend'
  }
  return(idx_sxt)
}

# --------------------------------------------------------------- #
# 3. Exposure layers at 25 km
# --------------------------------------------------------------- #
# Lookup of exposure files per agricultural dataset. Livestock layers are
# independent of the crop dataset. Fill/verify the paths before running.
exposure_files <- function(root, agr_dts){
  vln <- paste0(root,'/agex_results/agex_vulnerability')
  fls <- list(
    spam2010 = list(crop_diversity = paste0(vln,'/crop_classes_diversity_25km.tif'),
                    crop_vop       = paste0(vln,'/food_vop.tif')),
    spam2020 = list(crop_diversity = paste0(vln,'/crop_classes_diversity_25km_spam2020.tif'),
                    crop_vop       = paste0(vln,'/food_vop_spam2020.tif')) # CHECK: VoP 2020 file name
  )
  out <- fls[[agr_dts]]
  out$livestock_diversity <- paste0(vln,'/lsu_diversity_25km.tif')
  out$livestock_units     <- paste0(vln,'/lsu_total_25km.tif')
  return(out)
}

# Load the four exposure layers as one SpatRaster with standard names
load_exposures <- function(root, agr_dts){
  fls <- exposure_files(root, agr_dts)
  miss <- unlist(fls)[!file.exists(unlist(fls))]
  if(length(miss) > 0) stop(paste0('Missing exposure files:\n', paste(miss, collapse = '\n')))
  exp <- terra::rast(unlist(fls))
  names(exp) <- names(fls)
  return(exp)
}

# --------------------------------------------------------------- #
# 4. Hotspot identification (4x4 bivariate classification)
# --------------------------------------------------------------- #
# Cluster-level means of the extreme trend and each exposure layer.
cluster_metrics <- function(cls, idx_sxt, exp){
  names(cls) <- 'extreme_cluster'
  stk <- c(idx_sxt, exp)
  stk <- terra::resample(x = stk, y = cls, method = 'near') # ensure same grid
  mtr <- terra::zonal(x = stk, z = cls, fun = 'mean', na.rm = T)
  # Cluster area (km2) and number of pixels
  are <- terra::zonal(x = terra::cellSize(cls, unit = 'km') |> terra::mask(cls), z = cls, fun = 'sum', na.rm = T)
  names(are)[2] <- 'area_km2'
  npx <- terra::freq(cls); npx <- data.frame(extreme_cluster = npx$value, n_pixels = npx$count)
  mtr <- mtr |> dplyr::left_join(are, by = 'extreme_cluster') |> dplyr::left_join(npx, by = 'extreme_cluster')
  return(mtr)
}

# Quartile breaks as in 05 figure scripts: interior breaks from cluster-level
# means, lower/upper limits from pixel-level min/max.
quartile_breaks <- function(cl_values, px_values){
  brks <- stats::quantile(x = cl_values, probs = seq(0,1,1/4), na.rm = T)
  brks[1] <- min(px_values, na.rm = T)
  brks[5] <- max(px_values, na.rm = T)
  return(brks)
}

# Classify clusters into the 16 bivariate classes for every exposure variable.
# Returns list(classes = per-cluster table, thresholds = break table)
classify_hotspots <- function(mtr, idx_sxt, exp, hotspot_categories = c('4-4')){
  px_trd <- terra::values(idx_sxt, na.rm = T)
  trd_brks <- quartile_breaks(mtr$extreme_trend, px_trd)
  mtr$trend_class <- cut(x = mtr$extreme_trend, breaks = trd_brks, include.lowest = T) |> as.numeric()
  thr <- data.frame(variable = 'extreme_trend', q0 = trd_brks[1], q25 = trd_brks[2], q50 = trd_brks[3], q75 = trd_brks[4], q100 = trd_brks[5])
  for(v in names(exp)){
    px_v <- terra::values(exp[[v]], na.rm = T)
    v_brks <- quartile_breaks(mtr[[v]], px_v)
    mtr[[paste0(v,'_class')]] <- cut(x = mtr[[v]], breaks = v_brks, include.lowest = T) |> as.numeric()
    mtr[[paste0(v,'_bi_class')]] <- paste0(mtr$trend_class,'-',mtr[[paste0(v,'_class')]])
    mtr[[paste0(v,'_hotspot')]] <- mtr[[paste0(v,'_bi_class')]] %in% hotspot_categories
    thr <- rbind(thr, data.frame(variable = v, q0 = v_brks[1], q25 = v_brks[2], q50 = v_brks[3], q75 = v_brks[4], q100 = v_brks[5]))
  }
  # Composite hotspots as reported in the paper
  mtr$food_production_hotspot <- mtr$crop_diversity_hotspot | mtr$livestock_diversity_hotspot
  mtr$economic_value_hotspot  <- mtr$crop_vop_hotspot | mtr$livestock_units_hotspot
  rownames(thr) <- NULL
  return(list(classes = mtr, thresholds = thr))
}

# Rasterise a per-cluster logical/numeric column back to pixels
cluster_to_raster <- function(cls, tbl, column){
  names(cls) <- 'extreme_cluster'
  rcl <- cbind(tbl$extreme_cluster, as.numeric(tbl[[column]]))
  out <- terra::classify(x = cls, rcl = rcl, others = NA)
  names(out) <- column
  return(out)
}

# Full wrapper: one scenario = (index, agricultural dataset, clusters raster)
run_hotspots <- function(root, index, agr_dts, cls_file, scenario, out_dir){
  dir.create(path = out_dir, F, T)
  cls <- terra::rast(cls_file); names(cls) <- 'extreme_cluster'
  idx_sxt <- calc_extreme_trend(root, index)
  exp <- load_exposures(root, agr_dts)
  # Restrict pixel distributions to the pixels of this clusters' raster
  idx_sxt <- terra::resample(idx_sxt, cls, method = 'near') |> terra::mask(cls)
  exp     <- terra::resample(exp, cls, method = 'near') |> terra::mask(cls)
  mtr <- cluster_metrics(cls, idx_sxt, exp)
  res <- classify_hotspots(mtr, idx_sxt, exp)
  res$classes$scenario <- scenario; res$thresholds$scenario <- scenario
  utils::write.csv(x = res$classes, file = paste0(out_dir,'/agex_hotspots_',scenario,'.csv'), row.names = F)
  utils::write.csv(x = res$thresholds, file = paste0(out_dir,'/agex_hotspots_thresholds_',scenario,'.csv'), row.names = F)
  # Pixel-level layers used later for the comparisons
  hcols <- c('crop_diversity_hotspot','livestock_diversity_hotspot','crop_vop_hotspot','livestock_units_hotspot',
             'food_production_hotspot','economic_value_hotspot')
  hrs <- lapply(hcols, function(h) cluster_to_raster(cls, res$classes, h)) |> terra::rast()
  terra::writeRaster(x = c(cls, idx_sxt, exp, hrs), filename = paste0(out_dir,'/agex_hotspots_',scenario,'.tif'), overwrite = T)
  return(res)
}

# --------------------------------------------------------------- #
# 5. Comparison metrics between two scenarios
# --------------------------------------------------------------- #
# a and b: SpatRasters produced by run_hotspots (same grid).
# Returns one row per hotspot type.
compare_scenarios <- function(a, b, scenario_a, scenario_b, ctr = NULL,
                              hcols = c('crop_diversity_hotspot','livestock_diversity_hotspot',
                                        'crop_vop_hotspot','livestock_units_hotspot',
                                        'food_production_hotspot','economic_value_hotspot'),
                              scores = c('extreme_trend','crop_diversity','livestock_diversity','crop_vop','livestock_units'),
                              retain_thr = .5, seed = 1235, n_smp = 1e5){
  b <- terra::resample(b, a, method = 'near')
  are <- terra::cellSize(a[[1]], unit = 'km'); names(are) <- 'area'
  dfa <- terra::as.data.frame(c(are, a), cell = T, na.rm = F)
  dfb <- terra::as.data.frame(b, cell = T, na.rm = F)
  names(dfb)[-1] <- paste0(names(dfb)[-1],'_b')
  dfm <- dplyr::left_join(dfa, dfb, by = 'cell')
  cmn <- !is.na(dfm$extreme_cluster) & !is.na(dfm$extreme_cluster_b)

  # Partition similarity on common pixels
  ari <- aricode::ARI(dfm$extreme_cluster[cmn], dfm$extreme_cluster_b[cmn])
  nmi <- aricode::NMI(dfm$extreme_cluster[cmn], dfm$extreme_cluster_b[cmn])
  # Spearman correlation of cluster-level continuous scores at common pixels
  set.seed(seed)
  idx <- which(cmn); if(length(idx) > n_smp) idx <- sample(idx, n_smp)
  rho <- sapply(scores, function(s){
    stats::cor(dfm[idx, s], dfm[idx, paste0(s,'_b')], method = 'spearman', use = 'pairwise.complete.obs')
  })
  names(rho) <- paste0('spearman_',scores)

  out <- lapply(hcols, function(h){
    ha <- dfm[[h]] == 1 & !is.na(dfm[[h]])
    hb <- dfm[[paste0(h,'_b')]] == 1 & !is.na(dfm[[paste0(h,'_b')]])
    a_area <- sum(dfm$area[ha]); b_area <- sum(dfm$area[hb])
    i_area <- sum(dfm$area[ha & hb]); u_area <- sum(dfm$area[ha | hb])
    # Cluster-level retention: share of each baseline hotspot cluster's area
    # that is still hotspot in the alternative scenario
    ret <- data.frame(extreme_cluster = dfm$extreme_cluster[ha], area = dfm$area[ha], kept = hb[ha]) |>
      dplyr::group_by(extreme_cluster) |>
      dplyr::summarise(share = sum(area * kept) / sum(area), .groups = 'drop')
    data.frame(scenario_a = scenario_a, scenario_b = scenario_b, hotspot = h,
               n_clusters_a = length(unique(dfm$extreme_cluster[!is.na(dfm$extreme_cluster)])),
               n_clusters_b = length(unique(dfm$extreme_cluster_b[!is.na(dfm$extreme_cluster_b)])),
               common_pixels_pct = 100 * sum(cmn) / sum(!is.na(dfm$extreme_cluster)),
               ARI = ari, NMI = nmi,
               n_hotspots_a = length(unique(dfm$extreme_cluster[ha])),
               n_hotspots_b = length(unique(dfm$extreme_cluster_b[hb])),
               hotspot_area_a_km2 = a_area, hotspot_area_b_km2 = b_area,
               jaccard_area = ifelse(u_area > 0, i_area / u_area, NA),
               overlap_retained_pct = ifelse(a_area > 0, 100 * i_area / a_area, NA),
               clusters_retained_pct = ifelse(nrow(ret) > 0, 100 * mean(ret$share >= retain_thr), NA),
               t(rho))
  }) |> dplyr::bind_rows()

  # Regions where the classification changes (gained / lost hotspot area per country)
  if(!is.null(ctr)){
    cty <- terra::rasterize(x = ctr, y = a[[1]], field = 'adm0_iso')
    dfm$adm0_iso <- as.character(terra::as.data.frame(cty, na.rm = F)[,1]) # labels, not factor codes
    chg <- lapply(hcols, function(h){
      ha <- dfm[[h]] == 1 & !is.na(dfm[[h]])
      hb <- dfm[[paste0(h,'_b')]] == 1 & !is.na(dfm[[paste0(h,'_b')]])
      data.frame(adm0_iso = dfm$adm0_iso, area = dfm$area, lost = ha & !hb, gained = !ha & hb) |>
        dplyr::filter(lost | gained) |>
        dplyr::group_by(adm0_iso) |>
        dplyr::summarise(lost_km2 = sum(area[lost]), gained_km2 = sum(area[gained]), .groups = 'drop') |>
        dplyr::mutate(hotspot = h, scenario_a = scenario_a, scenario_b = scenario_b) |>
        dplyr::arrange(-(lost_km2 + gained_km2))
    }) |> dplyr::bind_rows()
    attr(out, 'regions') <- chg
  }
  return(out)
}
