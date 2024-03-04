## ------------------------------------------ ##
## Compute max-stable processes per cluster
## By: Harold Achicanoy
## WUR & ABC
## Dec. 2023
## ------------------------------------------ ##

options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
pacman::p_load(terra,clusterExtremes,landscapemetrics,geosphere,RColorBrewer,scales)

## Directories
root <- '//CATALOGUE/WFP_ClimateRiskPr1'

yrs <- 1980:2022

## Clustering results
# Cluster map
# fmado_r <- terra::rast(paste0(root,'/agroclimExtremes/results/clusters/ARG_fmado_clusters.tif'))
fmado_r <- terra::rast('D:/WORLD_fmado_trimmed_100_gs1.tif')
plot(fmado_r == 93)
fmado_r[fmado_r != 93] <- NA
fmado_r[!is.na(fmado_r)] <- 1; plot(fmado_r)
# fmado_r <- terra::focal(x = fmado_r, w = 3, fun = 'mean', na.rm = T); plot(fmado_r)

# Get coordinates
crd <- terra::as.data.frame(x = fmado_r, xy = T, cell = T, na.rm = T)

# Raw time series
cdd_roi <- terra::rast(paste0(root,'/agroclimExtremes/agex_indices/agex_cdd/agex_cdd_25km/one_s1_cdd_25km.tif'))
names(cdd_roi) <- paste0('Y',yrs)

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

# Get time series
tsr <- terra::extract(x = cdd_roi_fnl, y = crd[,c('x','y')], cell = T)

# Plotting time series
plot(x = 1979:2022, y = as.numeric(tsr[1,2:(ncol(tsr)-1)]), ty = 'l', col = scales::alpha('red', 0.05), ylim = c(-2,4))
for(i in 2:nrow(tsr)){
  lines(x = 1979:2022, y = as.numeric(tsr[i,2:(ncol(tsr)-1)]), col = scales::alpha('red', 0.05))
}; rm(i)

# Obtain medoid coordinates
med_crd <- t(apply(crd[,c('x','y')], 2, median)) |> base::as.data.frame()
plot(crd[,c('x','y')], pch = 20, col = 'black')
points(med_crd, pch = 20, col = 'red')
# Obtain medoid cell
med_cll <- terra::extract(fmado_r, med_crd, cell = T)$cell
# Get adjacent cells. TO DO. Think in a way to define the pixels to compute the Smith process
aux <- fmado_r
terra::values(aux) <- NA
aux[med_cll] <- 1
bfr <- terra::buffer(x = aux, width = 500000)
clls <- terra::extract(x = bfr, y = crd[,c('x','y')])
clls <- crd$cell[clls$layer]
# clls <- terra::adjacent(x = fmado_r, cells = med_cll, directions = '16') |> as.numeric()
# adjs <- terra::extract(x = fmado_r, y = terra::xyFromCell(object = fmado_r, cell = clls), cell = T) |> tidyr::drop_na() |> dplyr::pull(cell) |> as.numeric()

# Models' list
model_list <- lapply(1:10, function(i){
  
  # Sample 30 coordinates within the cluster of interest
  set.seed(i); smp_crd <- sample(x = clls, size = 30, replace = F)
  
  # Time series sampled data
  data_fit     <- t(as.matrix(tsr[tsr$cell %in% smp_crd,2:(ncol(tsr)-1)]))
  # Sampled coordinates
  coord_fit    <- as.matrix(crd[tsr$cell %in% smp_crd,c('x','y')])
  pair_weights <- get_pair_weights(data_fit, min_common_obs = 0)
  
  frech_bool <- F
  fitM <- tryCatch({
    if(frech_bool == T){
      fitM <- SpatialExtremes::fitmaxstab(data      = data_fit,
                                          coord     = coord_fit,
                                          marg.cov  = NULL,
                                          cov.mod   = 'gauss',
                                          iso       = F,
                                          weights   = pair_weights,
                                          fit.marge = F)
    } else {
      fitM <- SpatialExtremes::fitmaxstab(data       = data_fit,
                                          coord      = coord_fit,
                                          loc.form   <- loc ~ 1,
                                          scale.form <- scale ~ 1,
                                          shape.form <- shape ~ 1,
                                          marg.cov   = NULL,
                                          cov.mod    = 'gauss',
                                          iso        = F,
                                          weights    = pair_weights)
    }
    fitM;
  }, warning = function(w){
    cat("WARNING :",conditionMessage(w), "\n");
    return(NA)
  }, error = function(e){
    cat("ERROR :",conditionMessage(e), "\n");
    return(NA)
  })
  
  return(fitM)
  
})

utils_check_cov_ratio <- function(maxstable_model, level = 0.99){
  
  if(any(is.na(maxstable_model))){ return(NA) }
  
  param = maxstable_model$fitted.values
  cov11 = param[1]
  cov12 = param[2]
  cov22 = param[3]
  cov_mat = matrix(c(cov11, cov12, cov12, cov22), nrow = 2, byrow = TRUE)
  eigen_vals = eigen(cov_mat)$values
  r_values = 2*sqrt(eigen_vals*qchisq(level, 2))
  radius_ratio = r_values[2]/r_values[1]
  
  return(radius_ratio)
  
}

ratio_values <- lapply(model_list, utils_check_cov_ratio) |> unlist()

# This function gets the index i paramter from the max-stable model
utils_get_par_fun <- function(l, i){
  if(any(is.na(l))) return(NA)
  l$param[i]
}

# This function helps identify any suspect ellipses
utils_flag_ellipses <- function(maxstable_model_list, alpha = 0.05){
  
  # get parameters from list
  cov11 = lapply(maxstable_model_list, utils_get_par_fun, i = 1) |> unlist()
  cov12 = lapply(maxstable_model_list, utils_get_par_fun, i = 2) |> unlist()
  cov22 = lapply(maxstable_model_list, utils_get_par_fun, i = 3) |> unlist()
  param_info <- data.frame(sim_index = 1:length(maxstable_model_list),
                           cov11, cov12, cov22)
  
  # check ellipses
  check_ellipses <- param_info |>
    filter(cov11 < quantile(param_info$cov11, alpha/2, na.rm  = T) |
             cov11 > quantile(param_info$cov11, 1 - alpha/2, na.rm  = T) |
             cov22 < quantile(param_info$cov22, alpha/2, na.rm  = T) |
             cov22 > quantile(param_info$cov22, 1- alpha/2, na.rm  = T) |
             cov12 < quantile(param_info$cov12, alpha/2, na.rm  = T) |
             cov12 > quantile(param_info$cov12, 1- alpha/2, na.rm  = T))
  
  return(check_ellipses)
  
}
ellipse_alpha <- 0.05
# identify any suspect ellipses that do not fail ratio
check_ellipses <- utils_flag_ellipses(model_list, alpha = ellipse_alpha)
rerun_i <- check_ellipses$sim_index
model_list[rerun_i] <- NA

# get start values
start_list <- get_start_list(model_list, frech_bool)

convergence_issue = lapply(model_list, function(l){all(is.na(l))}) |> unlist() |> sum()
print(paste(convergence_issue, "samples with suspect convergence"))

iter <- 1
old_convergence_number <- convergence_issue

while(iter < max_iter){
  print("PARALLELISE ME")
  for(i in 1:num_samps){
    
    if(!any(is.na(model_list[[i]])) == T) next
    
    set.seed(i)
    set.seed(i); smp_crd <- sample(x = clls, size = 30, replace = F)
    
    # assume constant parameters across the space
    # less uncertainty that transforming to frechet
    # data_fit = frech_data[ , sample_stns] %>%
    #   as.matrix()
    # Time series sampled data
    data_fit     <- t(as.matrix(tsr[tsr$cell %in% smp_crd,2:(ncol(tsr)-1)]))
    # Sampled coordinates
    coord_fit    <- as.matrix(crd[tsr$cell %in% smp_crd,c('x','y')])
    pair_weights <- get_pair_weights(data_fit, min_common_obs = 0)
    
    fitM <- tryCatch({
      if(frech_bool == TRUE){
        fitM <- SpatialExtremes::fitmaxstab(data      = data_fit,
                                            coord     = coord_fit,
                                            marg.cov  = NULL,
                                            cov.mod   = 'gauss',
                                            iso       = F,
                                            weights   = pair_weights,
                                            fit.marge = F)
      }else{
        fitM <- SpatialExtremes::fitmaxstab(data       = data_fit,
                                            coord      = coord_fit,
                                            loc.form   <- loc ~ 1,
                                            scale.form <- scale ~ 1,
                                            shape.form <- shape ~ 1,
                                            marg.cov   = NULL,
                                            cov.mod    = "gauss",
                                            iso        = F,
                                            weights    = pair_weights)
      }
      fitM;
    }, warning = function(w){
      # cat("WARNING :",conditionMessage(e), "\n");
      return(NA)
    }, error = function(e){
      # cat("ERROR :",conditionMessage(e), "\n");
      return(NA)
    })
    
    model_list[[i]] = fitM
    
  }
  
  # -----------------------------------------------------------------------------
  
  # caluclate ratio of elliptical curves
  ratio_values <- lapply(model_list, utils_check_cov_ratio) |> unlist()
  
  rerun_i <- which(ratio_values < ratio_threshold)
  model_list[rerun_i] <- NA
  
  # -----------------------------------------------------------------------------
  
  # Don't want this step in the repeat iterations
  
  # # identify any suspect ellipses that do not fail ratio
  # check_ellipses <- utils_flag_ellipses(model_list, alpha = 0.1)
  # rerun_i = check_ellipses$sim_index
  # model_list[rerun_i] = NA
  
  # -----------------------------------------------------------------------------
  
  # get start values
  start_list <- get_start_list(model_list, frech_bool)
  
  convergence_issue <- lapply(model_list, function(l){all(is.na(l))}) |> unlist() |> sum()
  print(paste(convergence_issue, "samples with suspect convergence"))
  
  if(old_convergence_number == convergence_issue) iter = max_iter
  
  old_convergence_number = convergence_issue
  
  iter = iter + 1
  
}

# visualise the result
med_crd
ellipse_df <- get_ellipse_from_smith_model_list(model_list, medoid = med_crd)

# plot result
ggplot2::ggplot(data = ellipse_df) +
  ggplot2::geom_path(aes(x=x, y =y, group = sim_index), alpha = 0.25) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw()

crd$Values <- apply(tsr[,2:(ncol(tsr)-1)], 1, median)

crdref <- "+proj=longlat +datum=WGS84"

names(ellipse_df)[1:2] <- c('lon','lat')

ellipse_lst <- split(ellipse_df[,c('lon','lat')], ellipse_df$sim_index)
polygons    <- purrr::map(.x = ellipse_lst, .f = function(x){terra::vect(x = as.matrix(x), type = 'polygons', crs = crdref)})
polygons    <- terra::vect(polygons)
plot(polygons)
terra::writeVector(x = polygons, filename = 'D:/cluster_COL_test.gpkg', overwrite = T)

ggplot(data = crd) +
    geom_raster(aes(x = x, y = y, fill = Values)) +
    theme_void() +
    theme(legend.position = 'bottom') +
  coord_equal() +
  geom_path(data = ellipse_df, aes(x=x, y =y, group = sim_index), alpha = 0.25)

cdd_roi
chk <- fmado_r
chk[chk != 5] <- NA
chk[chk == 5] <- 1

cdd_roi_msk <- terra::mask(x = cdd_roi, mask = chk)
cdd_roi_msk[[1]]

tst <- cdd_roi_msk |> terra::as.data.frame(xy = T, cell = T)
tst$avg <- apply(tst[,paste0('Y',yrs)], 1, max)
tst$avg <- tst$Y1987
gg1 <- ggplot2::ggplot(data = tst, aes(x = x, y = y, colour = avg)) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::coord_equal()
ellipse_df$avg <- 1
gg1 +
  ggplot2::geom_path(data = ellipse_df, aes(x=x, y =y, group = sim_index), alpha = 1, size = 1.2, colour = 'red')
