## ------------------------------------------ ##
## Compute max-stable processes per cluster
## By: Harold Achicanoy
## WUR & ABC
## Dec. 2023
## ------------------------------------------ ##

options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
pacman::p_load(terra,clusterExtremes,landscapemetrics,geosphere,RColorBrewer,scales)

## Key arguments
root   <- '//CATALOGUE/WFP_ClimateRiskPr1'
yrs    <- 1980:2022
index  <- 'spei-6'
gs     <- 'one'
season <- 1

yrs <- 1980:2022

# Load extreme weather signatures
fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_results/clusters'), pattern = paste0('agex_global_',index,'_',gs,'_s',season,'_fmadogram_*.*.tif$'), full.names = T)
agex_sgn <- terra::rast(fls)
names(agex_sgn) <- 'extreme_signature'

sgntr <- 141

auxl_sgn <- agex_sgn
auxl_sgn[auxl_sgn != sgntr] <- NA
auxl_sgn[!is.na(auxl_sgn)] <- 1; plot(auxl_sgn)
# auxl_sgn <- terra::focal(x = auxl_sgn, w = 3, fun = 'mean', na.rm = F); plot(auxl_sgn)

# Get coordinates
crd <- terra::as.data.frame(x = auxl_sgn, xy = T, cell = T, na.rm = T) # Coordinates

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
tsr <- terra::extract(x = idx_roi_fnl, y = crd[,c('x','y')])
tsr$cell <- crd$cell

# Plotting time series
plot(x = 1980:2022, y = as.numeric(tsr[1,2:(ncol(tsr)-1)]), ty = 'l', col = scales::alpha('red', 0.05), ylim = c(-4,4))
for(i in 2:nrow(tsr)){
  lines(x = 1980:2022, y = as.numeric(tsr[i,2:(ncol(tsr)-1)]), col = scales::alpha('red', 0.05))
}; rm(i)

# Obtain medoid coordinates
med_crd <- t(apply(crd[,c('x','y')], 2, median)) |> base::as.data.frame()
plot(crd[,c('x','y')], pch = 20, col = 'black')
points(med_crd, pch = 20, col = 'red')
# Obtain medoid cell
med_cll <- terra::extract(auxl_sgn, med_crd, cell = T)$cell
# Get adjacent cells. TO DO. Think in a way to define the pixels to compute the Smith process
aux <- auxl_sgn
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
frech_bool <- F
start_list <- get_start_list(model_list, frech_bool)

convergence_issue = lapply(model_list, function(l){all(is.na(l))}) |> unlist() |> sum()
print(paste(convergence_issue, "samples with suspect convergence"))

frech_bool      = F
num_samps       = 10 # number of pixels to choose
samp_size       = 20 # number of samples per pixel
min_common_obs  = 0
min_pairs       = 0
ratio_threshold = 0.1
ellipse_alpha   = 0.05
max_iter        = 5

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

# Visualize the result
med_crd
ellipse_df <- get_ellipse_from_smith_model_list(model_list, medoid = med_crd)
ellipse_df <- ellipse_df[!is.na(ellipse_df$x),]

# plot result
ggplot2::ggplot(data = ellipse_df) +
  ggplot2::geom_path(aes(x=x, y =y, group = sim_index), alpha = 0.25) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw()

crd$Values <- apply(tsr[,2:(ncol(tsr)-1)], 1, max)

crdref <- "+proj=longlat +datum=WGS84"

names(ellipse_df)[1:2] <- c('lon','lat')

ellipse_lst <- split(ellipse_df[,c('lon','lat')], ellipse_df$sim_index)
polygons    <- purrr::map(.x = ellipse_lst, .f = function(x){terra::vect(x = as.matrix(x), type = 'polygons', crs = crdref)})
polygons    <- terra::vect(polygons)
plot(polygons)
terra::writeVector(x = polygons, filename = paste0(root,'/agroclimExtremes/agex_results/clusters/',index,'_',gs,'_s',season,'_signature_',sgntr,'.gpkg'), overwrite = T)

crd |>
  ggplot2::ggplot() +
  ggplot2::geom_raster(aes(x = x, y = y, fill = Values)) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = 'bottom') +
  ggplot2::coord_equal() +
  ggplot2::geom_path(data = ellipse_df, aes(x = lon, y = lat, group = sim_index), alpha = 0.25, size = 1.2)

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
