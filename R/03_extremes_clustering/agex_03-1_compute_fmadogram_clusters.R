## ------------------------------------------ ##
## Compute F-madogram distances
## By: Harold Achicanoy
## WUR & ABC
## Dec. 2023
## ------------------------------------------ ##

## R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra, geodata, Kendall, tidyverse, psych,
                                FactoMineR, factoextra, modifiedmk, NbClust,
                                rnaturalearth, RColorBrewer, MetBrewer,
                                fastcluster, eurostat, giscoR))

## Key arguments
root   <- '//CATALOGUE/WFP_ClimateRiskPr1'
yrs    <- 1980:2022
index  <- 'spei-6'
gs     <- 'two'
season <- 2

## List and load files
fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_25km'), pattern = paste0(gs,'_s',season,'_',index,'_'), full.names = T)
# Leave it open to compute for 10 km resolution
# fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_10km'), pattern = paste0(gs,'_s',season,'_',index,'_'), full.names = T)
idx <- terra::rast(fls)
names(idx) <- paste0('Y',yrs)

# Mask by resampled growing season id
ngs <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_nseasons_25km.tif'))
if(gs == 'one'){ngs[ngs == 2] <- NA} else {ngs[ngs == 1] <- NA; ngs[!is.na(ngs)] <- 1}

## Region of interest
## Crop by region of interest
idx_roi <- terra::mask(x = idx, mask = ngs)
# # Country level
# roi <- geodata::gadm(country = 'AUS', level = 0, path = tempdir(), version = 'latest')
# # Continental level
# roi <- terra::vect(paste0(root,'/agroclimExtremes/agex_raw_data/agex_africa_shape.gpkg'))
# roi <- eurostat::get_eurostat_geospatial(output_class = 'sf', resolution = '60', year = '2016', crs = '4326')
# roi <- terra::vect(roi)
# idx_roi <- terra::crop(x = idx, y = terra::ext(roi))
# idx_roi <- terra::mask(x = idx_roi, mask = roi)

# Maybe it is important to show a map of trends, positive and negative significant trends

## Final index: applying transformations
transformations <- function(x){
  if(!all(is.na(x))){
    ## Remove NA and -Inf values
    # x[is.na(x)] <- 0
    # x[is.infinite(x)] <- 0
    # Scale the time series
    x_scl <- scale(x)
    # Mann-Kendall test for evaluating trend significance
    trd <- Kendall::MannKendall(x_scl)$sl
    cnd <- ifelse(trd <= 0.05, T, F)
    if(cnd){ # Remove trend by fitting loess regression
      dfm <- data.frame(x = 1:length(x_scl), y = x_scl)
      fit <- loess(y ~ x, data = dfm)
      x_ntrd <- dfm$y - fit$fitted
    } else { # Same scaled time series
      x_ntrd <- x_scl
    }
  } else {
    x_ntrd <- rep(NA, length(x))
  }
  return(x_ntrd)
}

idx_roi_fnl <- terra::app(x = idx_roi, fun = function(i, ff) ff(i), cores = 20, ff = transformations)
names(idx_roi_fnl) <- paste0('Y',yrs); gc(T)

## Transform raster into a data.frame of stationary processes
idx_roi_ntrd <- terra::as.data.frame(x = idx_roi_fnl, xy = T, cell = T, na.rm = T)
x <- t(idx_roi_ntrd[,4:ncol(idx_roi_ntrd)])

## Compute F-madogram distance
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
  # DD_fmado = dist(t(V), method = "manhattan", diag = TRUE, upper = TRUE)/(2*Tnb)
  DD_fmado <- parallelDist::parDist(x = t(V), method = 'manhattan', diag = T, upper = T)/(2*Tnb)
  return(DD_fmado)
}
# Trim distribution of F-madogram distances greater than 1/6
cap_fmado_dist <- function(DD_fmado){
  DD_fmado_adj <- pmin(DD_fmado, rep(1/6, length(DD_fmado)))
  if(any(DD_fmado_adj > 1/6, na.rm = TRUE) == 1)
    stop("ERROR: cap_fmado_dist didn't work")
  return(DD_fmado_adj)
}

## Obtain distance matrices
fmado_dist <- get_fmado_dist(x); gc(T)          # F-madogram distances
fmado_dist <- cap_fmado_dist(fmado_dist); gc(T) # Truncated F-madogram distances
# saveRDS(object = fmado_dist, file = paste0(root,'/agroclimExtremes/agex_results/WORLD_',index,'_',gs,'_s',season,'_fmado_distance.rds'))
eucld_dist <- parallelDist::parDist(x = t(x),   # Euclidean distances matrix
                                    method = 'euclidean',
                                    diag = T, upper = T); gc(T)

# Hierarchical clustering for distance matrices
fmado_hcl <- fastcluster::hclust(fmado_dist, method = 'ward.D2') # F-madogram
eucld_hcl <- fastcluster::hclust(eucld_dist, method = 'ward.D2') # Euclidean
gc(T)

# Query data.frame
qry <- idx_roi_ntrd[,c('cell','x','y')]

# Optimal number of clusters
# 1. Execute a PAM (Partition Around Medoids) algorithm with a sample of pixels
# 2. Vary the number of clusters between 2 to 200 and test each partition
# 3. Evaluate a metric (e.g., silhouette) and pick the optimal K-number
source('https://raw.githubusercontent.com/haachicanoy/agroclimExtremes/main/R/03_extremes_clustering/agex_03-2_optimal_number_clusters.R')

# 4. Execute the hierarchical clustering with that K-number
qry$fmado_cluster <- cutree(tree = fmado_hcl, k = optimal_k)
qry$eucld_cluster <- cutree(tree = eucld_hcl, k = optimal_k)

# Raster template
tmp <- idx_roi_fnl[[1]]
terra::values(tmp) <- NA

# Assign clusters to pixels
fmado_r <- tmp; fmado_r[qry$cell] <- qry$fmado_cluster; plot(fmado_r)
eucld_r <- tmp; eucld_r[qry$cell] <- qry$eucld_cluster; plot(eucld_r)

terra::writeRaster(x = fmado_r, filename = paste0(root,'/agroclimExtremes/agex_results/clusters/agex_global_',index,'_',gs,'_s',season,'_fmadogram_k',optimal_k,'.tif'), overwrite = T)
terra::writeRaster(x = eucld_r, filename = paste0(root,'/agroclimExtremes/agex_results/clusters/agex_global_',index,'_',gs,'_s',season,'_euclidean_k',optimal_k,'.tif'), overwrite = T)
