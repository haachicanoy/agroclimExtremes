# -------------------------------------------------------- #
# Transnational potential mapping (supplementary material)
# Political and agro-ecological zones comparison
# Crops and livestock clusters
# By: Harold Achicanoy
# WUR & ABC
# Oct 2024
# -------------------------------------------------------- #

options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
pacman::p_load(terra, tidyverse)

root <- '//CATALOGUE/AgroclimExtremes'

# Extreme drought clusters
shp <- terra::vect(paste0(root,'/agex_results/agex_results_clusters/vct_agex_global_spei-6.gpkg'))
shp <- shp[,'extreme_cluster']

## Drought clusters span political borders ----
wrl <- terra::vect('D:/worldshape.gpkg') # Country's shapefile
wrl <- wrl[wrl$name_en != 'Antarctica',]
wrl_dsv <- terra::aggregate(wrl, dissolve = T) # Remove country borders
wrl_lns <- terra::as.lines(wrl)         # From polygons to lines
wrl_dsv_lns <- terra::as.lines(wrl_dsv) # From polygons to lines
tst <- terra::erase(x = wrl_lns, y = wrl_dsv_lns); plot(tst) # Remove external lines from the shapefile
shp_brd <- terra::relate(shp, tst, 'intersects') # Drought clusters intersected with country borders
trnt_cls <- shp[apply(shp_brd, 1, sum) > 0,] # Transnational clusters
plot(trnt_cls)
terra::writeVector(x = trnt_cls, filename = paste0(root,'/agex_results/agex_results_clusters/vct_agex_transnational.gpkg'), overwrite = T)
nt_cls <- shp[apply(shp_brd, 1, sum) == 0,] # National clusters
plot(nt_cls)
terra::writeVector(x = nt_cls, filename = paste0(root,'/agex_results/agex_results_clusters/vct_agex_national.gpkg'), overwrite = T)

## Hotspots for crops ----
crop_hotspots <- utils::read.csv(paste0(root,'/agex_results/agex_crops_hotspots.csv'))
crop_shp <- terra::merge(x = shp, y = crop_hotspots)
terra::writeVector(x = crop_shp, filename = paste0(root,'/agex_results/agex_results_clusters/vct_agex_crops_hotspots.gpkg'), overwrite = T)

## Hotspots for livestock ----
lvst_hotspots <- utils::read.csv(paste0(root,'/agex_results/agex_livestock_hotspots.csv'))
lvst_shp <- terra::merge(x = shp, y = lvst_hotspots)
terra::writeVector(x = lvst_shp, filename = paste0(root,'/agex_results/agex_results_clusters/vct_agex_livestock_hotspots.gpkg'), overwrite = T)

## Drought clusters span agro-ecological borders ----
agex_sgn_num <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn_num) <- 'extreme_cluster'

gaez <- terra::rast(paste0(root,'/agex_raw_data/gaez_v4/aez/aez_v9v2_CRUTS32_Hist_8110_100_avg.tif'))
gaez <- gaez |> terra::crop(terra::ext(agex_sgn_num)) |> terra::resample(agex_sgn_num, method = 'near')
gaez_smp <- c(2, 4, 37) # 6, 26, 31
gaez_tls <- c('Tropics, lowland',
              'Tropics, lowland',
              # 'Tropics, lowland',
              # 'Sub-tropics, cool',
              # 'Temperate, moderate',
              'Temperate, cool')
gaez_sbt <- c('Semi-arid, with soil/terrain limitations',
              'Sub-humid, with soil/terrain limitations',
              # 'Humid, with soil/terrain limitations',
              # 'Semi-arid, with soil/terrain limitations',
              # 'Dry, no soil/terrain limitations',
              'Dry, no soil/terrain limitations')
gaez_ltr <- c('a) ','b) ','c) ')

# Color palette definition
n_col    <- length(unique(terra::values(agex_sgn_num,na.rm = T))) # Number of unique colors needed
col_pltt <- MetBrewer::met.brewer(name = 'Signac', n = n_col) # Color palette
set.seed(1235); col_pltt <- sample(x = col_pltt, size = n_col, replace = F); rm(n_col) # Randomize colors assignment

# Coordinates
agex_sgn_dfm <- terra::as.data.frame(x = agex_sgn_num, xy = T, cell = T)
agex_sgn_dfm$extreme_cluster <- as.factor(agex_sgn_dfm$extreme_cluster)
# Setting-up colors
stp <- data.frame(extreme_cluster = factor(sort(unique(terra::values(agex_sgn_num, na.rm = T)))), color = col_pltt)
agex_sgn_dfm <- dplyr::left_join(x = agex_sgn_dfm, y = stp, by = 'extreme_cluster')
agex_sgn_dfm$extreme_cluster <- as.factor(agex_sgn_dfm$extreme_cluster)
sub_stp <- unique(agex_sgn_dfm[,c('extreme_cluster','color')])

gg <- list()
for(i in 1:length(gaez_smp)){
  gaez_tmp <- (gaez == gaez_smp[i])
  gaez_tmp_shp <- terra::as.polygons(x = gaez_tmp); names(gaez_tmp_shp) <- 'Filter'
  gaez_tmp_shp <- gaez_tmp_shp[gaez_tmp_shp$Filter == 1,]
  aux <- terra::mask(x = agex_sgn_num, mask = gaez_tmp_shp)
  n_ext <- terra::ext(terra::trim(aux))
  n_ext$xmin <- -180; n_ext$xmax <- 180
  aux_dfm <- terra::as.data.frame(x = aux, xy = T, cell = T)
  wrl_dsv_aux <- terra::crop(x = wrl_dsv, y = n_ext)
  
  ggm <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = agex_sgn_dfm[which(agex_sgn_dfm$cell %in% aux_dfm$cell),], aes(x, y, fill = extreme_cluster)) +
    ggplot2::geom_sf(data = sf::st_as_sf(wrl_dsv_aux), fill = NA) +
    ggplot2::coord_sf() +
    ggplot2::scale_fill_manual(values = sub_stp$color, breaks = sub_stp$extreme_cluster) +
    ggplot2::theme_classic() +
    ggplot2::xlab('') +
    ggplot2::ylab('') +
    ggplot2::ggtitle(label = paste0(gaez_ltr[i], gaez_tls[i]), subtitle = gaez_sbt[i]) +
    ggplot2::theme(legend.position = 'none',
                   text = element_text(family = 'serif'),
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank())
  gg[[i]] <- ggm
}

layout <- matrix(1:3, nrow = 3, ncol = 1, byrow = T)
fig <- gridExtra::grid.arrange(gg[[1]], gg[[2]], gg[[3]], layout_matrix = layout)
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/supplementary/FigureS_agex_spanning_gaez.png'), plot = fig, device = 'png', width = 5, height = 6, units = 'in', dpi = 350)
