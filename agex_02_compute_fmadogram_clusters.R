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
fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_indices/agex_cdd'), pattern = 'one_s[0-9]_cdd_', full.names = T)
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
                              if(cnd){ # Remove trend
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
names(cdd_roi_fnl) <- paste0('Y',yrs)

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

## Obtain distance matrices
geogr <- as.matrix(scale(cdd_roi_ntrd[,c('x','y')]))
geogr_dist <- parallelDist::parDist(x = geogr, # Geographical distance
                                    method = 'euclidean',
                                    diag = T, upper = T)
fmado_dist  <- get_fmado_dist(x)       # F-madogram distances matrix
fmado_dist2 <- fmado_dist + geogr_dist # F-madogram + geography
fmado_dist2 <- .6 * fmado_dist + .4 * geogr_dist # F-madogram + geography
eucld_dist  <- parallelDist::parDist(x = t(x), # Euclidean distances matrix
                                     method = 'euclidean',
                                     diag = T, upper = T)
eucld_dist2 <- eucld_dist + geogr_dist # Euclidean + geography

# Hierarchical clustering for distance matrices
# fmado_hcl <- hclust(fmado_dist, method = 'complete') # F-madogram
# eucld_hcl <- hclust(eucld_dist, method = 'complete') # Euclidean
fmado_hdb <- dbscan::hdbscan(x = fmado_dist, minPts = 100)
fmado_hcl <- fastcluster::hclust(fmado_dist, method = 'complete') # F-madogram
# single, complete, average, mcquitty, ward.D, ward.D2, centroid, median
fmado_hcl_sng <- fastcluster::hclust(fmado_dist, method = 'single') # F-madogram
fmado_hcl_cmp <- fastcluster::hclust(fmado_dist, method = 'complete') # F-madogram
fmado_hcl_avg <- fastcluster::hclust(fmado_dist, method = 'average') # F-madogram
fmado_hcl_mcq <- fastcluster::hclust(fmado_dist, method = 'mcquitty') # F-madogram
fmado_hcl_wrd <- fastcluster::hclust(fmado_dist, method = 'ward.D') # F-madogram
fmado_hcl_wd2 <- fastcluster::hclust(fmado_dist, method = 'ward.D2') # F-madogram
fmado_hcl_cnt <- fastcluster::hclust(fmado_dist, method = 'centroid') # F-madogram
fmado_hcl_mdn <- fastcluster::hclust(fmado_dist, method = 'median') # F-madogram
eucld_hcl <- fastcluster::hclust(eucld_dist, method = 'complete') # Euclidean
geogr_hcl <- fastcluster::hclust(geogr_dist, method = 'complete') # F-madogram

# Query data.frame
qry <- cdd_roi_ntrd[,c('cell','x','y')]

qry$fmado_cluster <- cutree(tree = fmado_hcl, k = 30)
qry$fmado_cluster_sng <- cutree(tree = fmado_hcl_sng, k = 30)
qry$fmado_cluster_cmp <- cutree(tree = fmado_hcl_cmp, k = 30)
qry$fmado_cluster_avg <- cutree(tree = fmado_hcl_avg, k = 30)
qry$fmado_cluster_mcq <- cutree(tree = fmado_hcl_mcq, k = 30)
qry$fmado_cluster_wrd <- cutree(tree = fmado_hcl_wrd, k = 30)
qry$fmado_cluster_wd2 <- cutree(tree = fmado_hcl_wd2, k = 30)
qry$fmado_cluster_cnt <- cutree(tree = fmado_hcl_cnt, k = 30)
qry$fmado_cluster_mdn <- cutree(tree = fmado_hcl_mdn, k = 30)
qry$eucld_cluster <- cutree(tree = eucld_hcl, k = 30)

saveRDS(object = qry, file = paste0(root,'/agroclimExtremes/results/cluster_info.rds'))

res <- NbClust::NbClust(diss     = fmado_dist,
                        distance = NULL,
                        min.nc   = 2,
                        max.nc   = 140,
                        method   = 'ward.D2',
                        index    = 'frey')

# cut_heights = data.frame(h = full_tree$height,
#                          k = (nrow(cdd_roi_dfm)-1):1)

# Raster template
tmp <- cdd_roi_fnl[[1]]
terra::values(tmp) <- NA

fmado_r <- tmp; fmado_r[qry$cell] <- qry$fmado_cluster; plot(fmado_r)
fmado_r_sng <- tmp; fmado_r_sng[qry$cell] <- qry$fmado_cluster_sng
fmado_r_cmp <- tmp; fmado_r_cmp[qry$cell] <- qry$fmado_cluster_cmp
fmado_r_avg <- tmp; fmado_r_avg[qry$cell] <- qry$fmado_cluster_avg
fmado_r_mcq <- tmp; fmado_r_mcq[qry$cell] <- qry$fmado_cluster_mcq
fmado_r_wrd <- tmp; fmado_r_wrd[qry$cell] <- qry$fmado_cluster_wrd
fmado_r_wd2 <- tmp; fmado_r_wd2[qry$cell] <- qry$fmado_cluster_wd2
fmado_r_cnt <- tmp; fmado_r_cnt[qry$cell] <- qry$fmado_cluster_cnt
fmado_r_mdn <- tmp; fmado_r_mdn[qry$cell] <- qry$fmado_cluster_mdn
eucld_r <- tmp; eucld_r[qry$cell] <- qry$eucld_cluster; plot(eucld_r)

fmado_r <- c(fmado_r_cmp,
             fmado_r_avg,
             fmado_r_mcq,
             fmado_r_wrd,
             fmado_r_wd2)

fmado_r <- terra::rast('//CATALOGUE/WFP_ClimateRiskPr1/agroclimExtremes/agex_results/clusters/test.tif')
names(fmado_r) <- c('complete','average','mcq','ward','ward2')

pca_model <- terra::prcomp(x = fmado_r)
pca_r <- terra::predict(fmado_r, pca_model)

terra::writeRaster(x = fmado_r, filename = paste0(root,'/agroclimExtremes/agex_results/clusters/EUR_fmado_clusters.tif'), overwrite = T)
terra::writeRaster(x = eucld_r, filename = paste0(root,'/agroclimExtremes/agex_results/clusters/EUR_eucld_clusters.tif'), overwrite = T)
saveRDS(object = cdd_roi_ntrd, file = paste0(root,'/agroclimExtremes/results/ts_processed.rds'))

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
