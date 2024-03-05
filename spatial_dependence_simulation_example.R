library(pacman)
pacman::p_load(clusterExtremes, scales, dbscan)

# Generate an arbitary set of 600 points (stations)
set.seed(1)
n.site    <- 600
locations <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(locations) <- c('lon','lat')
dim(locations)
plot(locations)

# Create a test clustering and classification dataframe
# make two disjoint regions (two clusters)
# make some of the probabilities less than 0.5
knn_info <- data.frame(x = locations[,1], y = locations[,2]) |>
  dplyr::mutate(knn_id = 1) |>
  dplyr::mutate(cluster_id = 1) |>
  dplyr::mutate(prob = 1) |>
  dplyr::mutate(knn_id = if_else((x < 6 & y < 6) | (x > 7 & y > 7), 2, knn_id)) |>
  dplyr::mutate(cluster_id = if_else((x < 6 & y < 6) | (x > 7 & y > 7), 2, cluster_id)) |>
  dplyr::mutate(prob = if_else(x < 5.5 & y < 5.5, prob, 0.4))
dim(knn_info)
head(knn_info)
plot(knn_info$x, knn_info$y, col = knn_info$knn_id, pch = 20)

# Generate some max-stable process data (time series data for stations)
set.seed(2)
data <- SpatialExtremes::rmaxstab(n       = 100,
                                  coord   = locations,
                                  cov.mod = 'gauss',
                                  cov11   = 1,
                                  cov12   = 0,
                                  cov22   = 1)
dim(data)
plot(as.numeric(data[,1]), ty = 'l', col = scales::alpha('forestgreen', .05), ylim = c(0, 850))
for(i in 2:ncol(data)){
  lines(as.numeric(data[,i]), col = scales::alpha('forestgreen',.05))
}; rm(i)

# Create grid covering the coordinates (don't work)
# grid <- get_grid_for_classification(coords = knn_info %>% select(x,y), grid_space = 0.25)
# so, running step by step
grid_space <- 0.25
min_dist   <- 1.5 * grid_space
full_grid  <- generate_grid(coords = knn_info |> select(x,y), grid_space = grid_space)
near_grid  <- utils_reduce_grid(coords = knn_info |> select(x,y), full_grid = full_grid, min_dist = min_dist)
plot(full_grid); points(near_grid, col = 2, pch = 20); points(locations, col = 4, pch = 20)
grid <- near_grid
names(grid) <- c('x','y')
dim(grid)
# so grid object is the near grid that is gonna be used. regular points that are
# close to the stations data

# Mimic the classification of the knn_info to the grid object
grid <- grid |>
  dplyr::mutate(knn_id = 1) |>
  dplyr::mutate(knn_id = if_else((x < 6 & y < 6) | (x > 7 & y > 7), 2, knn_id))
dim(grid)
plot(grid$x, grid$y, col = grid$knn_id, pch = 20)

# fit the smith model (don't work)
# model_list <- fit_smith_model(data = data, knn_info = knn_info,
#                               medoid_indexes = c(100, 601), use_id = 2,
#                               grid, min_class_prob = 0.5, high_class_prob = 0.75,
#                               frech_bool = TRUE,
#                               num_samps = 10, samp_size = 20,
#                               min_common_obs = 0, min_pairs = 0,
#                               ratio_threshold = 0.1, ellipse_alpha = 0.05,
#                               max_iter = 3,
#                               save_output = FALSE, output_dir = NULL,
#                               reference_id = NULL)
# medoid_indexes = c(100, 601)
use_id          = 1    # Cluster 1?
min_class_prob  = 0.5  # Something related with classification
high_class_prob = 0.75 # Something related with classification
frech_bool      = T    # Parameter of max-stable model
num_samps       = 10   # number of pixels to choose
samp_size       = 20   # number of samples per pixel
min_common_obs  = 0
min_pairs       = 0
ratio_threshold = 0.1
ellipse_alpha   = 0.05
max_iter        = 3
save_output     = F

fit_smith_model

# Checks
id_check = all(knn_info$cluster_id %in% knn_info$knn_id) &
  all(knn_info$knn_id %in% knn_info$cluster_id)
if(id_check == FALSE) warning("Error: cluster_id and knn_id does not match")

id1_check = all(knn_info$knn_id %in% grid$knn_id) &
  all(grid$knn_id %in% knn_info$knn_id)
if(id1_check == FALSE) stop("Error: grid_id and knn_id does not match")

range_check = all(max(grid$x) >= knn_info$x) &
  all(min(grid$x) <= knn_info$x) &
  all(max(grid$y) >= knn_info$y) &
  all(min(grid$y) <= knn_info$y)
if(range_check == FALSE) warning("Warning: grid does not cover range of station data")

if(is.integer(num_samps)) stop("Error: num_samps must be an integer")
if(is.integer(samp_size)) stop("Error:samp_size must be an integer")
if(is.integer(max_iter)) stop("Error: max_iter must be an integer")

if(min_class_prob > 1 | min_class_prob < 0)
  stop("Error: incorrect value for min_class_prob specified")

if(ellipse_alpha > 0.5 | ellipse_alpha < 0)
  stop("Error: incorrect value for ellipse_alpha specified")

if(ncol(data) != nrow(knn_info))
  stop("Error: There data and knn_info have incompatible dimensions")

# create a new class_id that is binary
fit_info <- knn_info |> dplyr::mutate(class_id = (knn_id == use_id))
grid <- grid |> dplyr::mutate(class_id = (knn_id == use_id))
plot(fit_info[fit_info$class_id,c('x','y')]); points(grid[grid$class_id,c('x','y')], col = 2, pch = 20)

# check the region is connected
# connection_check <- check_clusters_connected(coords = fit_info |> select(x, y, class_id), grid = grid)
# doesn't work
coords <- fit_info |> dplyr::select(x, y, class_id)
grid
grid_space       <- min(dist(grid |> dplyr::select(x, y)))
all_coords       <- rbind(coords |> dplyr::select(x, y), grid |> dplyr::select(x,y))
class_coords     <- rbind(coords |> dplyr::select(class_id), grid |> dplyr::select(class_id))
dd               <- dist(all_coords)
dd_class         <- dist(cbind(class_coords, class_coords))
dd_sum           <- dd + dd_class
db               <- dbscan::dbscan(dd_sum, eps = grid_space, minPts = 1)
db_id            <- db$cluster[1:nrow(coords)]
connection_check <- db_id

# get possible stations for fitting
possible_stns <- which(connection_check == 2 &
                         fit_info$prob > min_class_prob &
                         !((fit_info$cluster_id != fit_info$knn_id) & fit_info$prob > high_class_prob)) |> as.numeric()

if(length(possible_stns) < samp_size + 1){
  warning("Warning: too few stations for sampling, reduce samp_size")
  samp_size = length(possible_stns) - 1
}

model_list <- vector('list',num_samps)
for(i in 1:num_samps){
  
  set.seed(i)
  sample_stns = sample(possible_stns, samp_size, replace = F)
  
  data_fit = data[,sample_stns] |> as.matrix()
  
  coord_fit = fit_info[sample_stns,] |> dplyr::select(x,y) |> as.matrix()
  
  pair_weights = get_pair_weights(data_fit, min_common_obs)
  if(sum(pair_weights) < min_pairs){
    warning("Warning: too few common pairs for fitting, reduce sample size")
  }
  
  fitM <- tryCatch({
    if(frech_bool == TRUE){
      fitM <- fitmaxstab(data      = data_fit,
                         coord     = coord_fit,
                         marg.cov  = NULL,
                         cov.mod   = 'gauss',
                         iso       = F,
                         weights   = pair_weights,
                         fit.marge = F)
    }else{
      fitM <- fitmaxstab(data       = data_fit,
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
  
  model_list[[i]] <- fitM
  # fitM
  
}

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

# calculate ratio of elliptical curves
ratio_values <- lapply(model_list, utils_check_cov_ratio) |> unlist()

rerun_i <- which(ratio_values < ratio_threshold)
model_list[rerun_i] <- NA

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

# identify any suspect ellipses that do not fail ratio
check_ellipses <- utils_flag_ellipses(model_list, alpha = ellipse_alpha)
rerun_i = check_ellipses$sim_index
model_list[rerun_i] = NA

# get start values
start_list = get_start_list(model_list, frech_bool)

convergence_issue = lapply(model_list, function(l){all(is.na(l))}) |> unlist() |> sum()
print(paste(convergence_issue, "samples with suspect convergence"))

iter = 1
old_convergence_number = convergence_issue

while(iter < max_iter){
  print("PARALLELISE ME")
  for(i in 1:num_samps){
    
    if(!any(is.na(model_list[[i]])) == T) next
    
    set.seed(i)
    sample_stns <- sample(possible_stns, samp_size, replace = F)
    
    # assume constant parameters across the space
    # less uncertainty that transforming to frechet
    # data_fit = frech_data[ , sample_stns] %>%
    #   as.matrix()
    data_fit <- data[,sample_stns] |> as.matrix()
    
    coord_fit <- fit_info[sample_stns,] |> select(x,y) |> as.matrix()
    
    pair_weights <- get_pair_weights(data_fit, min_common_obs)
    
    fitM <- tryCatch({
      if(frech_bool == TRUE){
        fitM <- fitmaxstab(data      = data_fit,
                           coord     = coord_fit,
                           marg.cov  = NULL,
                           cov.mod   = 'gauss',
                           iso       = F,
                           weights   = pair_weights,
                           fit.marge = F)
      }else{
        fitM <- fitmaxstab(data       = data_fit,
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

model_list

# visualise the result
medoid_coords = data.frame(x = 3, y = 3)
ellipse_df <- get_ellipse_from_smith_model_list(model_list, medoid = medoid_coords)

# plot result
ell_plot <- ggplot(data = ellipse_df) +
  geom_path(aes(x = x, y = y, group = sim_index), alpha = 0.25) +
  coord_equal() +
  theme_bw()
ell_plot
