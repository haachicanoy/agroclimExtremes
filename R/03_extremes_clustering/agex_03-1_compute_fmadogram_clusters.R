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

## Directories
root <- '//CATALOGUE/WFP_ClimateRiskPr1'

yrs <- 1979:2022

## List and load files
fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_indices/agex_cdd/agex_cdd_10km'), pattern = 'one_s[0-9]_cdd_', full.names = T)
cdd <- terra::rast(fls)
names(cdd) <- paste0('Y',yrs)

## Region of interest
# roi <- geodata::gadm(country = 'AUS', level = 0, path = tempdir(), version = 'latest')
# roi <- terra::vect(paste0(root,'/agroclimExtremes/agex_raw_data/agex_africa_shape.gpkg'))
roi <- eurostat::get_eurostat_geospatial(output_class = 'sf', resolution = '60', year = '2016', crs = '4326')
roi <- terra::vect(roi)

# ## Graph of one index-year
# wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
# wrl <- terra::aggregate(wrl)
# # my.palette <- RColorBrewer::brewer.pal(n = 20, name = 'Set1')
# my.palette <- MetBrewer::met.brewer(name = 'Tam', n = 20)
# tmp <- cdd[[terra::nlyr(cdd)]]
# tmp <- terra::trim(x = tmp)
# hist(terra::values(tmp))
# qnt <- quantile(x = terra::values(tmp), probs = 0.98, na.rm = T)
# terra::values(tmp)[terra::values(tmp) > qnt] <- qnt
# png(filename = "D:/test.png", width = 3132, height = 2359, units = 'px')
# plot(wrl, ext = terra::ext(tmp))
# plot(tmp, add = T, col = my.palette, plg = list(cex = 5), cex.main = 7)
# dev.off()

## Crop by region of interest
cdd_roi <- terra::crop(x = cdd, y = terra::ext(roi))
cdd_roi <- terra::mask(x = cdd_roi, mask = roi)

## Apply transformations
cdd_roi_fnl <- terra::app(x = cdd_roi,
                          fun = function(x){
                            if(!all(is.na(x))){
                              # Remove NA and -Inf values
                              x[is.na(x)] <- 0
                              x[is.infinite(x)] <- 0
                              # Scale the time series
                              x_scl <- scale(x)
                              # Mann-Kendall test for evaluating trend significance
                              trd <- Kendall::MannKendall(x_scl)$sl
                              cnd <- ifelse(trd <= 0.05, T, F)
                              if(cnd){ # Remove trend by fitting loess regression
                                dfm <- data.frame(x = yrs, y = x_scl)
                                fit <- loess(y ~ x, data = dfm)
                                x_ntrd <- dfm$y - fit$fitted
                              } else { # Same scaled time series
                                x_ntrd <- x_scl
                              }
                            } else {
                              x_ntrd <- rep(NA, length(x))
                            }
                            return(x_ntrd)
                          })
names(cdd_roi_fnl) <- paste0('Y',yrs); gc(T)

## Transform into a data.frame
cdd_roi_ntrd <- terra::as.data.frame(x = cdd_roi_fnl, xy = T, cell = T, na.rm = T)
x <- t(cdd_roi_ntrd[,4:ncol(cdd_roi_ntrd)])

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
fmado_dist <- get_fmado_dist(x); gc(T) # F-madogram distances
# Saving F-madogram distance matrix. Not efficient.
# saveRDS(object = fmado_dist, file = paste0(root,'/agroclimExtremes/agex_results/cdd_extreme_distance.rds'))
fmado_dist <- cap_fmado_dist(fmado_dist); gc(T) # Truncated F-madogram distances
eucld_dist  <- parallelDist::parDist(x = t(x), # Euclidean distances matrix
                                     method = 'euclidean',
                                     diag = T, upper = T); gc(T)

# Hierarchical clustering for distance matrices
fmado_hcl <- fastcluster::hclust(fmado_dist, method = 'ward.D2') # F-madogram
eucld_hcl <- fastcluster::hclust(eucld_dist, method = 'ward.D2') # Euclidean
# fmado_hdb <- dbscan::hdbscan(x = fmado_dist, minPts = 100) # Hierarchical DBSCAN
# geogr_hcl <- fastcluster::hclust(geogr_dist, method = 'complete') # F-madogram

# Query data.frame
qry <- cdd_roi_ntrd[,c('cell','x','y')]

# How to determine number of clusters?
# 1. Execute a K-means/PAM algorithm with a sample of pixels
# 2. Vary the number of clusters between 2 to 200/100 to 200
# 3. Evaluate a metric (e.g., silhouette) and pick a K-number
# 4. Execute the hierarchical clustering with that K-number

qry$fmado_cluster <- cutree(tree = fmado_hcl, k = 10)
qry$eucld_cluster <- cutree(tree = eucld_hcl, k = 10)

# saveRDS(object = qry, file = paste0(root,'/agroclimExtremes/results/cluster_info.rds'))

# res <- NbClust::NbClust(diss     = fmado_dist,
#                         distance = NULL,
#                         min.nc   = 2,
#                         max.nc   = 140,
#                         method   = 'ward.D2',
#                         index    = 'frey')

# cut_heights = data.frame(h = full_tree$height,
#                          k = (nrow(cdd_roi_dfm)-1):1)

# Raster template
tmp <- cdd_roi_fnl[[1]]
terra::values(tmp) <- NA

fmado_r <- tmp; fmado_r[qry$cell] <- qry$fmado_cluster; plot(fmado_r)
eucld_r <- tmp; eucld_r[qry$cell] <- qry$eucld_cluster; plot(eucld_r)

# fmado_r <- c(fmado_r_cmp,
#              fmado_r_avg,
#              fmado_r_mcq,
#              fmado_r_wrd,
#              fmado_r_wd2)
# 
# fmado_r <- terra::rast('//CATALOGUE/WFP_ClimateRiskPr1/agroclimExtremes/agex_results/clusters/test.tif')
# names(fmado_r) <- c('complete','average','mcq','ward','ward2')
# 
# pca_model <- terra::prcomp(x = fmado_r)
# pca_r <- terra::predict(fmado_r, pca_model)

terra::writeRaster(x = fmado_r, filename = paste0(root,'/agroclimExtremes/agex_results/clusters/AUS_cdd_fmado_trimmed_clusters.tif'), overwrite = T)
terra::writeRaster(x = eucld_r, filename = paste0(root,'/agroclimExtremes/agex_results/clusters/AUS_cdd_eucld_clusters.tif'), overwrite = T)
# saveRDS(object = cdd_roi_ntrd, file = paste0(root,'/agroclimExtremes/results/ts_processed.rds'))

par(mfrow = c(2,5))
for(i in 1:10){
  plot(roi)
  plot(terra::mask(x = fmado_r, mask = fmado_r, maskvalues = i, inverse = T), add = T)
}; rm(i)

# Describe results
cdd_sts <- cdd_roi_ntrd[,c('cell','x','y')]
cdd_sts$avg <- apply(cdd_roi_ntrd[,4:ncol(cdd_roi_ntrd)], 1, mean)
cdd_sts$std <- apply(cdd_roi_ntrd[,4:ncol(cdd_roi_ntrd)], 1, sd)
cdd_sts$mdn <- apply(cdd_roi_ntrd[,4:ncol(cdd_roi_ntrd)], 1, median)
cdd_sts$min <- apply(cdd_roi_ntrd[,4:ncol(cdd_roi_ntrd)], 1, min)
cdd_sts$max <- apply(cdd_roi_ntrd[,4:ncol(cdd_roi_ntrd)], 1, max)
cdd_sts$trd <- apply(cdd_roi_ntrd[,4:ncol(cdd_roi_ntrd)], 1, function(x) trend::sens.slope(x)$estimates)

psych::describeBy(x = cdd_sts, group = qry$fmado_cluster)
pca.res <- FactoMineR::PCA(X = cdd_sts[,4:ncol(cdd_sts)], scale.unit = T, graph = F)
factoextra::fviz_pca_ind(pca.res,  label = 'none', habillage = qry$fmado_cluster)
