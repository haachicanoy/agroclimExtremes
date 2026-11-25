# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Sensitivity analysis: Adjusted Rand Index
# By: Harold Achicanoy
# WUR & ABC
# Created in August 2026
# Modified in September 2026
# --------------------------------------------------------------- #

# R options and user-defined functions ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,tidyverse,mclust))
list.files2 <- Vectorize(FUN = list.files, vectorize.args = 'path')
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')

root   <- '//CATALOGUE/AgroclimExtremes'


# Extreme clusters for places with one growing seasons
agex_sgn <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_one_s1_fmadogram_clean.tif'))
# Lower limit: minimum number of pixels per extreme cluster based on gs1
lwr_lmt <- min(table(as.numeric(terra::values(agex_sgn, na.rm = T))))

agex_sgn_s1 <- terra::rast("//CATALOGUE/AgroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-3_two_s1_fmadogram_k222.tif")
agex_sgn_s2 <- terra::rast("//CATALOGUE/AgroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-3_two_s2_fmadogram_k200.tif")

# Counting the number of coincident extreme clusters
dfm <- terra::as.data.frame(x = c(agex_sgn_s1, agex_sgn_s2), xy = T, cell = T)
names(dfm)[(ncol(dfm)-1):ncol(dfm)] <- c('s1','s2')
dfm$combined <- paste0(dfm$s1,'-',dfm$s2)
table(dfm$combined)

# Remove clusters with a fewer number of pixels
aux <- base::as.data.frame(sort(table(dfm$combined), decreasing = T))
names(aux) <- c('combined','count')
aux$extreme_cluster <- 1:nrow(aux)
aux$extreme_cluster[aux$count < lwr_lmt] <- NA

dfm <- dplyr::left_join(x = dfm, y = aux, by = 'combined')

agex_sgn_aux <- agex_sgn; rm(agex_sgn)
terra::values(agex_sgn_aux) <- NA
agex_sgn_aux[dfm$cell] <- dfm$extreme_cluster

agex_sgn_cln <- agex_sgn_aux
n_clusters   <- length(unique(as.numeric(terra::values(agex_sgn_aux, na.rm = T))))
for(cl in 1:n_clusters){
  
  cat(paste0('Processing CLUSTER: ',cl,'\n'))
  aux <- agex_sgn_aux
  aux[aux != cl] <- NA
  
  # Get cluster's coordinates
  aux_dfm <- terra::as.data.frame(x = aux, xy = T, cell = T, na.rm = T)
  
  # Get outliers from coordinates
  set.seed(1235)
  out_pre <- OutliersO3::O3prep(data = aux_dfm[,c('x','y')], method = c('BAC','adjOut','DDC','MCD')) # 'PCS'
  out_msg <- capture.output(try(OutliersO3::O3plotM(out_pre), silent = T))
  if(length(out_msg) > 1){
    out_res <- OutliersO3::O3plotM(out_pre)
    
    # Number of outliers per method
    out_cns <- out_res$nOut
    print(out_cns)
    n_out <- getmode(out_cns)
    
    cat(paste0('Outliers identified: ',n_out,'\n\n'))
    if(n_out > 0){
      out_cnd <- out_res$outsTable$Method %in% names(which(out_res$nOut == n_out))
      out_css <- unique(out_res$outsTable[out_cnd,'Case'])
      agex_sgn_cln[aux_dfm[out_css,'cell']] <- NA
    }
  }
  
}

outfile <- paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-3_two_fmadogram_clean.tif')
terra::writeRaster(x = agex_sgn_cln, outfile, overwrite = T)


library(mclust)

compare_ari <- function(reference, alternative) {
  valid <- complete.cases(reference, alternative)
  
  adjustedRandIndex(
    reference[valid],
    alternative[valid]
  )
}

compare_ari(reference = spei6, alternative = spei3)

tst <- terra::xapp(x = spei6, y = spei3, fun = compare_ari)

pairwise_membership_agreement <- function(
    cluster_a,
    cluster_b,
    weights = NULL,
    n_pairs = 1e6,
    seed = 123
) {
  
  valid <- which(
    !is.na(cluster_a) &
      !is.na(cluster_b)
  )
  
  if (is.null(weights)) {
    weights <- rep(1, length(cluster_a))
  }
  
  prob <- weights[valid]
  prob <- prob / sum(prob)
  
  set.seed(seed)
  
  i <- sample(valid, n_pairs, replace = TRUE, prob = prob)
  j <- sample(valid, n_pairs, replace = TRUE, prob = prob)
  
  same_a <- cluster_a[i] == cluster_a[j]
  same_b <- cluster_b[i] == cluster_b[j]
  
  mean(same_a == same_b)
}

tst <- pairwise_membership_agreement(cluster_a = spei6, cluster_b = spei3)

cluster_overlap_matrix <- function(
    cluster_ref,
    cluster_alt,
    cell_weights
) {
  
  valid <- complete.cases(
    cluster_ref,
    cluster_alt,
    cell_weights
  )
  
  a <- cluster_ref[valid]
  b <- cluster_alt[valid]
  w <- cell_weights[valid]
  
  ref_ids <- sort(unique(a))
  alt_ids <- sort(unique(b))
  
  out <- matrix(
    NA_real_,
    nrow = length(ref_ids),
    ncol = length(alt_ids),
    dimnames = list(ref_ids, alt_ids)
  )
  
  area_ref <- tapply(w, a, sum)
  area_alt <- tapply(w, b, sum)
  
  for (i in seq_along(ref_ids)) {
    for (j in seq_along(alt_ids)) {
      
      intersection <- sum(
        w[a == ref_ids[i] & b == alt_ids[j]]
      )
      
      union <- area_ref[i] +
        area_alt[j] -
        intersection
      
      out[i, j] <- intersection / union
    }
  }
  
  out
}

tst <- cluster_overlap_matrix(cluster_ref = spei6, cluster_alt = spei3)

# r <- terra::rast("//CATALOGUE/AgroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_one_s1_fmadogram_clean.tif")

terra::rast("//CATALOGUE/AgroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif")
terra::rast("//CATALOGUE/AgroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-3_combined_fmadogram_clean.tif")
terra::rast("//CATALOGUE/AgroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-12_combined_fmadogram_clean.tif")

r <- terra::rast("//CATALOGUE/AgroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-3_combined_fmadogram_clean.tif")
s <- terra::rast("//CATALOGUE/AgroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-12_combined_fmadogram_clean.tif")

dfm1 <- terra::as.data.frame(c(r,s), xy = T, cell = T)
dfm2 <- terra::as.data.frame(c(r,s), xy = T, cell = T, na.rm = T)

names(dfm1)[c(4:5)] <- c("SPEI6", "SPEI3")
names(dfm2)[c(4:5)] <- c("SPEI6", "SPEI3")

mclust::adjustedRandIndex(dfm1$SPEI6, dfm1$SPEI3)
mclust::adjustedRandIndex(dfm2$SPEI6, dfm2$SPEI3)
