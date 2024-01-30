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
                                FactoMineR, factoextra, modifiedmk, NbClust))

## Directories
root <- '//CATALOGUE/WFP_ClimateRiskPr1'

yrs <- 1981:2019

## List and load files
fls <- list.files(path = paste0(root,'/agroclimExtremes/indices/cdd'), pattern = 'one_s[0-9]_cdd_', full.names = T)
cdd <- terra::rast(fls)
names(cdd) <- paste0('Y',yrs)

## Region of interest
roi <- geodata::gadm(country = 'ARG', level = 0, path = tempdir(), version = 'latest')
# roi <- terra::vect(paste0(root,'/agroclimExtremes/data/shapefiles/africa.gpkg'))

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
# F-madogram distances matrix
fmado_dist <- get_fmado_dist(x)
# Euclidean distances matrix
eucld_dist <- parallelDist::parDist(x = t(x), method = 'euclidean', diag = T, upper = T)

# Hierarchical clustering for distance matrices
fmado_hcl <- hclust(fmado_dist, method = 'complete') # F-madogram
eucld_hcl <- hclust(eucld_dist, method = 'complete') # Euclidean

# Query data.frame
qry <- cdd_roi_ntrd[,c('cell','x','y')]

qry$fmado_cluster <- cutree(tree = fmado_hcl, k = 10)
qry$eucld_cluster <- cutree(tree = eucld_hcl, k = 10)

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
eucld_r <- tmp; eucld_r[qry$cell] <- qry$eucld_cluster; plot(eucld_r)

terra::writeRaster(x = fmado_r, filename = paste0(root,'/agroclimExtremes/results/clusters/ARG_fmado_clusters.tif'), overwrite = T)
terra::writeRaster(x = fmado_r, filename = paste0(root,'/agroclimExtremes/results/clusters/ARG_eucld_clusters.tif'), overwrite = T)
saveRDS(object = cdd_roi_ntrd, file = paste0(root,'/agroclimExtremes/results/ts_processed.rds'))

par(mfrow = c(2,5))
for(i in 1:10){
  plot(roi)
  plot(terra::mask(x = fmado_r, mask = fmado_r, maskvalues = i, inverse = T), add = T)
}; rm(i)

# Describe results
psych::describeBy(x = cdd_roi_ntrd[,paste0('Y',yrs)], group = qry$fmado_cluster)
pca.res <- FactoMineR::PCA(X = cdd_per_dfm[,paste0('Y',1981:2019)], scale.unit = T, ncp = length(1981:2019), graph = F)
factoextra::fviz_pca_ind(pca.res,  label = 'none', habillage = cdd_per_dfm$cluster)
