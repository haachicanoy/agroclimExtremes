options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra, geodata, Kendall, tidyverse, psych,
                                FactoMineR, factoextra, modifiedmk, NbClust,
                                rnaturalearth, RColorBrewer, MetBrewer,
                                fastcluster, eurostat, giscoR))

yrs <- 1980:2022
idx_roi <- terra::rast('/home/ruser/data/one_s1_spei-6_10km.tif')
names(idx_roi) <- paste0('Y',yrs)

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
# Query data.frame
qry <- idx_roi_ntrd[,c('cell','x','y')]

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

## Obtain distance matrices
fmado_dist <- get_fmado_dist(x); gc(T)          # F-madogram distances
# fmado_dist <- cap_fmado_dist(fmado_dist); gc(T) # Truncated F-madogram distances
rm(x); gc(T)

# Optimal number of clusters
# 1. Execute a PAM (Partition Around Medoids) algorithm with a sample of pixels
# 2. Vary the number of clusters between 2 to 200 and test each partition
# 3. Evaluate a metric (e.g., silhouette) and pick the optimal K-number
source('https://raw.githubusercontent.com/haachicanoy/agroclimExtremes/main/R/03_extremes_clustering/agex_03-2_optimal_number_clusters.R')
rm(idx_roi_ntrd, idx_roi_fnl); gc(T)

optimal_k

# Hierarchical clustering for distance matrices
fmado_hcl <- fastcluster::hclust(fmado_dist, method = 'ward.D2'); gc(T) # F-madogram

# 4. Execute the hierarchical clustering with that K-number
qry$fmado_cluster <- cutree(tree = fmado_hcl, k = optimal_k)

# Raster template
tmp <- idx_roi_fnl[[1]]
terra::values(tmp) <- NA

# Assign clusters to pixels
fmado_r <- tmp; fmado_r[qry$cell] <- qry$fmado_cluster; plot(fmado_r)

index  <- 'spei-6'
gs     <- 'one'
season <- 1

terra::writeRaster(x = fmado_r, filename = paste0('/home/ruser/data/agex_global_',index,'_',gs,'_s',season,'_fmadogram_k',optimal_k,'.tif'), overwrite = T)

saveRDS(fmado_dist, '/home/ruser/data/fmado_dist.rds')

fmado_dist_mtrx <- as.matrix(fmado_dist)
saveRDS(fmado_dist_mtrx, '/home/ruser/data/fmado_dist_mtrx.rds')
