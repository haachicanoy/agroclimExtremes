## ------------------------------------------ ##
## Find optimal number of clusters
## By: Harold Achicanoy
## WUR & ABC
## Feb. 2024
## ------------------------------------------ ##

## R options and packages loading
suppressMessages(pacman::p_load(cluster,future,furrr))

# Sample a number of pixels
set.seed(1245)
prc <- ifelse(gs == 'one', 0.01, 0.1)
smp <- sample(x = idx_roi_ntrd$cell, size = prc * dim(idx_roi_ntrd)[1], replace = F)
# Sampled locations
x_smp <- t(idx_roi_ntrd[idx_roi_ntrd$cell %in% smp,4:ncol(idx_roi_ntrd)])

# Obtain distance matrices
fmado_smp_dist <- get_fmado_dist(x_smp); gc(T)   # F-madogram distances
fmado_smp_dist <- cap_fmado_dist(fmado_smp_dist) # Truncated F-madogram distances
gc(T)

future::plan(cluster, workers = 40, gc = T)
sil_width <- furrr::future_map(.x = 2:200, .f = function(k){
  model <- cluster::pam(x = fmado_smp_dist, k = k)
  results <- data.frame(k = k, avg_silhouette = model$silinfo$avg.width)
  return(results)
}) |> dplyr::bind_rows()
future:::ClusterRegistry('stop'); gc(T)

# # Plot the relationship between k and avg_silhouette
# ggplot(sil_width, aes(x = k, y = avg_silhouette)) +
#   geom_line() + geom_point() +
#   scale_x_continuous(breaks = 2:200)

optimal_k <- sil_width$k[which.max(sil_width$avg_silhouette)]
rm(fmado_smp_dist, x_smp, smp, sil_width, prc)
