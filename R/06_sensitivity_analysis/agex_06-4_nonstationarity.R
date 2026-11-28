# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Sensitivity analysis: non-stationarity of the extremal dependence
# structure (Reviewer 3.1)
# By: Harold Achicanoy
# WUR & ABC
# Created in November 2026
# --------------------------------------------------------------- #
# Two checks per growing-season setting (one s1, two s1, two s2):
# A. Raw vs detrended: cluster the scaled raw series (no Loess detrending)
#    with the same K and compare it with the detrended partition (ARI).
#    Answers whether cluster boundaries from residuals diverge from those
#    of the raw data.
# B. Split-half with a block-permutation null: cluster 1980-2000 and
#    2001-2022 separately and compute the ARI between both partitions.
#    Under a stationary dependence structure, a chronological split should
#    be as similar as a random split of the same years. The null distribution
#    is built by splitting the years at random in contiguous blocks
#    (preserving short-range serial dependence). Empirical p-value =
#    share of random splits with ARI <= chronological ARI.
#    Both chronological and random splits use the same full-period detrended
#    series, so the only difference is WHICH years enter each half.
# Note: with ~21 values per half, F-madogram estimates are noisy; this is
# exactly why the observed ARI must be read against the null, not alone.
# --------------------------------------------------------------- #

# R options and user-defined functions ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra, dplyr, aricode))
source('https://raw.githubusercontent.com/haachicanoy/agroclimExtremes/refs/heads/main/R/06_sensitivity_analysis/agex_06-0_sensitivity_functions.R')

# Key arguments ----
root     <- '//CATALOGUE/AgroclimExtremes'
out_dir  <- paste0(root,'/agex_results/agex_sensitivity/nonstationarity')
index    <- 'spei-6'
yrs      <- 1980:2022
brk_yr   <- 2000     # first half: 1980-2000 (21 yrs); second half: 2001-2022 (22 yrs)
n_boot   <- 1000    # random block splits for the null distribution
blk_len  <- 3        # block length (years)
smp_prc  <- 0.10     # share of pixels used in check B (same sample in all runs)
dir.create(out_dir, F, T)

# Random split of the years in contiguous blocks
block_split <- function(n_yrs, n_first, blk_len){
  blk <- split(1:n_yrs, ceiling(seq_along(1:n_yrs)/blk_len))
  ord <- unname(unlist(blk[sample(length(blk))]))
  first <- sort(ord[1:n_first])
  list(first = first, second = setdiff(1:n_yrs, first))
}

stp <- data.frame(gs = c('one','two','two'), season = c(1,1,2))

res <- lapply(1:nrow(stp), function(i){

  gs <- stp$gs[i]; season <- stp$season[i]
  cat(paste0('Processing: ',gs,' s',season,'\n'))
  k <- get_baseline_k(root, index, gs, season)

  idx_roi <- load_index_ts(root, index, gs, season, yrs)
  # Detrended (as in 03-1) and raw (scaled only) series
  # idx_dtr <- terra::app(x = idx_roi, fun = function(i, ff) ff(i), cores = 20, ff = function(x) transformations(x, detrend = TRUE))
  # idx_raw <- terra::app(x = idx_roi, fun = function(i, ff) ff(i), cores = 20, ff = function(x) transformations(x, detrend = FALSE))
  idx_dtr <- terra::app(x = idx_roi, fun = function(i, tf) {library(Kendall); tf(i, detrend = TRUE)}, cores = 20, tf = transformations)
  idx_raw <- terra::app(x = idx_roi, fun = function(i, tf) {library(Kendall); tf(i, detrend = FALSE)}, cores = 20, tf = transformations)
  dfm_dtr <- terra::as.data.frame(x = idx_dtr, cell = T, na.rm = T)
  dfm_raw <- terra::as.data.frame(x = idx_raw, cell = T, na.rm = T)
  dfm_raw <- dfm_raw[match(dfm_dtr$cell, dfm_raw$cell),]
  x_dtr <- t(dfm_dtr[,-1]); x_raw <- t(dfm_raw[,-1]); gc(T)

  # A. Raw vs detrended (all pixels, same K) ----
  cls_dtr <- fmado_cluster(x_dtr, k)
  cls_raw <- fmado_cluster(x_raw, k)
  ari_raw <- aricode::ARI(cls_dtr, cls_raw)
  r <- idx_roi[[1]]; terra::values(r) <- NA
  r_raw <- r; r_raw[dfm_dtr$cell] <- cls_raw; names(r_raw) <- 'extreme_cluster'
  terra::writeRaster(x = r_raw, filename = paste0(out_dir,'/agex_global_',index,'_',gs,'_s',season,'_fmadogram_k',k,'_raw.tif'), overwrite = T)
  # Share of pixels with a significant trend (context for interpreting A)
  sig <- apply(x_raw, 2, function(z) Kendall::MannKendall(z)$sl <= 0.05)

  # B. Split-half vs block-permutation null (pixel sample) ----
  set.seed(1235)
  smp <- sort(sample(1:ncol(x_dtr), size = round(smp_prc * ncol(x_dtr))))
  xs <- x_dtr[,smp]
  n_first <- sum(yrs <= brk_yr)
  h1 <- fmado_cluster(xs[1:n_first,], k)
  h2 <- fmado_cluster(xs[(n_first+1):length(yrs),], k)
  ari_chr <- aricode::ARI(h1, h2)
  ari_h1_full <- aricode::ARI(h1, cls_dtr[smp])
  ari_h2_full <- aricode::ARI(h2, cls_dtr[smp])

  set.seed(2026)
  ari_null <- sapply(1:n_boot, function(b){
    sp <- block_split(length(yrs), n_first, blk_len)
    aricode::ARI(fmado_cluster(xs[sp$first,], k), fmado_cluster(xs[sp$second,], k))
  })
  utils::write.csv(x = data.frame(gs = gs, season = season, boot = 1:n_boot, ARI = ari_null),
                   file = paste0(out_dir,'/agex_nonstationarity_null_',index,'_',gs,'_s',season,'.csv'), row.names = F)

  data.frame(index = index, gs = gs, season = season, k = k,
             n_pixels = ncol(x_dtr), pct_pixels_sig_trend = 100 * mean(sig, na.rm = T),
             ARI_raw_vs_detrended = ari_raw,
             n_pixels_split = length(smp),
             ARI_chronological_halves = ari_chr,
             ARI_first_half_vs_full = ari_h1_full,
             ARI_second_half_vs_full = ari_h2_full,
             ARI_null_median = stats::median(ari_null),
             ARI_null_q05 = stats::quantile(ari_null, .05),
             ARI_null_q95 = stats::quantile(ari_null, .95),
             p_value = mean(ari_null <= ari_chr))
}) |> dplyr::bind_rows()

res[,-(1:3)] <- round(res[,-(1:3)], 3)
utils::write.csv(x = res, file = paste0(out_dir,'/agex_TableS_nonstationarity_',index,'.csv'), row.names = F)
res

# Combined raw partition (same merge rules as 03-4/03-5 can be applied with
# the rasters written above if a map is needed for the Supplementary material).
