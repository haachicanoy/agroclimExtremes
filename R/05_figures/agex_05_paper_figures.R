## ------------------------------------------ ##
## Paper figures/maps
## By: Harold Achicanoy
## WUR & ABC
## May 2024
## ------------------------------------------ ##

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,MetBrewer,rnaturalearth,tidyverse,ggpubr,tseries,biscale,pals,cowplot,gridExtra,grid))

# Activate Windows fonts
extrafont::font_import()
extrafont::loadfonts(device = 'win')
extrafont::fonts()

# list.files2 <- Vectorize(FUN = list.files, vectorize.args = 'path')
# grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')
# get_trend <- function(x){
#   if(!all(is.na(x))){
#     x <- as.numeric(na.omit(x))
#     y <- as.numeric(trend::sens.slope(x)$estimates)
#   } else { y <- NA }
#   return(y)
# }
# get_ext_trend <- function(x){
#   if(!all(is.na(x))){
#     x <- as.numeric(na.omit(x))
#     dfm <- data.frame(time = 1:length(x), ts = x)
#     qrf <- quantreg::rq(formula = ts ~ time, tau = .90, data = dfm)
#     y <- as.numeric(qrf$coefficients[2])
#   } else { y <- NA }
#   return(y)
# }
# getmode <- function(v){
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }
# cv_fun <- function(x, na.rm = T){
#   sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm) * 100
# }

## Key arguments
root   <- '//CATALOGUE/WFP_ClimateRiskPr1'
index  <- 'spei-6'

### Extreme drought clusters (categorical)
agex_sgn <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn) <- 'extreme_cluster'
cls <- data.frame(extreme_cluster = sort(unique(terra::values(agex_sgn))))
cls$id <- cls$extreme_cluster
cls <- cls[,c('id','extreme_cluster')]
cls$extreme_cluster <- as.character(cls$extreme_cluster)
levels(agex_sgn) <- cls; rm(cls)

### Extreme drought clusters (numerical)
agex_sgn_num <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn_num) <- 'extreme_cluster'

### Extreme drought clusters metrics
agex_sgn_metrics <- utils::read.csv(paste0(root,'/agroclimExtremes/agex_results/agex_all_metrics.csv'))

### Load global shapefile
wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
wrl <- wrl[,c('adm0_iso','name_en','continent','subregion')]

## Shapefiles by continent. Some of them cropped
afr <- rnaturalearth::ne_countries(scale = 'large', continent = 'africa', returnclass = 'sv')
eur <- rnaturalearth::ne_countries(scale = 'large', continent = 'europe', returnclass = 'sv')
eur <- terra::crop(x = eur, y = terra::ext(c(-20,100,30,75))) # xmin, xmax = 150, ymin, ymax
asi <- rnaturalearth::ne_countries(scale = 'large', continent = 'asia', returnclass = 'sv')
oce <- rnaturalearth::ne_countries(scale = 'large', continent = 'oceania', returnclass = 'sv')
oce <- terra::crop(x = oce, y = terra::ext(c(90,180,-54.75,20.55)))
nam <- rnaturalearth::ne_countries(scale = 'large', continent = 'north america', returnclass = 'sv')
nam <- terra::crop(x = nam, y = terra::ext(c(-140,-50,-0.39,60)))
sam <- rnaturalearth::ne_countries(scale = 'large', continent = 'south america', returnclass = 'sv')

# Put them all together
shp <- list(nam, sam, afr, eur, asi, oce); rm(nam, sam, afr, eur, asi, oce)

# # For masking the rasters
# agex_sgn_afr <- terra::crop(x = agex_sgn, y = terra::ext(afr)); # agex_sgn_afr <- terra::mask(x = agex_sgn_afr, mask = afr)
# agex_sgn_eur <- terra::crop(x = agex_sgn, y = terra::ext(eur)); # agex_sgn_eur <- terra::mask(x = agex_sgn_eur, mask = eur)
# agex_sgn_asi <- terra::crop(x = agex_sgn, y = terra::ext(asi)); # agex_sgn_asi <- terra::mask(x = agex_sgn_asi, mask = asi)
# agex_sgn_oce <- terra::crop(x = agex_sgn, y = terra::ext(oce)); # agex_sgn_oce <- terra::mask(x = agex_sgn_oce, mask = oce)
# agex_sgn_nam <- terra::crop(x = agex_sgn, y = terra::ext(nam)); # agex_sgn_nam <- terra::mask(x = agex_sgn_nam, mask = nam)
# agex_sgn_sam <- terra::crop(x = agex_sgn, y = terra::ext(sam)); # agex_sgn_sam <- terra::mask(x = agex_sgn_sam, mask = sam)

# To test
# myTheme <- rasterVis::rasterTheme(region = col_pltt)
# # Europe
# rasterVis::levelplot(x            = agex_sgn_eur,
#                      layers       = 'extreme_cluster',
#                      par.settings = myTheme,
#                      colorkey     = F,
#                      maxpixels    = terra::ncell(agex_sgn))

## Figure 1
## Extreme drought clusters per continent
# Color palette definition
n_col    <- length(unique(terra::values(agex_sgn,na.rm = T))) # Number of unique colors needed
col_pltt <- MetBrewer::met.brewer(name = 'Signac', n = n_col) # Color palette
set.seed(1); col_pltt <- sample(x = col_pltt, size = n_col, replace = F); rm(n_col) # Randomize colors assignment

# Create a list of maps (per continent)
ggm <- list()
for(i in 1:length(shp)){
  
  # Coordinates
  agex_sgn_dfm <- terra::as.data.frame(x = terra::crop(agex_sgn, terra::ext(shp[[i]])), xy = T, cell = T)
  # Setting-up colors
  stp <- data.frame(extreme_cluster = factor(sort(unique(terra::values(agex_sgn, na.rm = T)))), color = col_pltt)
  agex_sgn_dfm <- dplyr::left_join(x = agex_sgn_dfm, y = stp, by = 'extreme_cluster')
  agex_sgn_dfm$extreme_cluster <- as.factor(agex_sgn_dfm$extreme_cluster)
  sub_stp <- unique(agex_sgn_dfm[,c('extreme_cluster','color')])
  # Map
  ggm[[i]] <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = agex_sgn_dfm, aes(x, y, fill = extreme_cluster)) +
    ggplot2::geom_sf(data = sf::st_as_sf(shp[[i]]), fill = NA) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = 'none', legend.key.width = unit(1.3, 'cm')) +
    ggplot2::scale_fill_manual(values = sub_stp$color, breaks = sub_stp$extreme_cluster)
  
}

# Arrange them into one graph and add panel labels -----
fig1 <- ggpubr::ggarrange(ggm[[1]], ggm[[2]], ggm[[3]], ggm[[4]], ggm[[5]], ggm[[6]],
                             labels = c('a)', 'b)', 'c)', 'd)', 'e)', 'f)'),
                             ncol = 2, nrow = 3, font.label = list(size = 30, family = 'serif', face = 'plain'))
fig1
ggplot2::ggsave(filename = paste0('D:/Figure1_paper1.png'), plot = fig1, device = 'png', width = 10, height = 12.5, units = 'in', dpi = 350); rm(fig1)

## Figure 2 ------------------------
## Percentage of extreme drought clusters in one, two, or more countries
### Collaboration potential -----------------------
col_pltt <- MetBrewer::met.brewer('Ingres', n = 3)
fig2 <- agex_sgn_metrics |>
  dplyr::group_by(countries_count, growing_seasons) |>
  dplyr::summarize(n()) |>
  dplyr::rename(Count = `n()`) |>
  dplyr::ungroup() |>
  dplyr::group_by(growing_seasons) |>
  dplyr::mutate(Freq = Count/sum(Count) * 100) |>
  dplyr::ungroup() |>
  dplyr::mutate(growing_seasons = factor(growing_seasons)) |>
  ggplot2::ggplot(aes(x = countries_count, y = Freq, fill = growing_seasons)) +
  ggplot2::geom_bar(stat = 'identity', position = position_dodge()) +
  ggplot2::scale_x_continuous(breaks = 1:13) +
  ggplot2::scale_fill_manual(values = c(col_pltt[1], col_pltt[2])) +
  ggplot2::xlab('Number of countries') +
  ggplot2::ylab('Percentage of extreme drought clusters (%)') +
  ggplot2::labs(fill = 'Number of growing\nseasons') +
  ggplot2::theme_minimal() +
  ggplot2::theme(text            = element_text(size = 17, colour = 'black', family = 'serif'),
                 axis.text.x     = element_text(size = 16, colour = 'black'),
                 axis.text.y     = element_text(size = 16, colour = 'black'),
                 axis.title      = element_text(size = 20, colour = 'black'),
                 legend.text     = element_text(size = 13, colour = 'black'),
                 legend.title    = element_text(size = 15, colour = 'black'),
                 plot.title      = element_text(size = 25, colour = 'black'),
                 plot.subtitle   = element_text(size = 17, colour = 'black'),
                 strip.text.x    = element_text(size = 17, colour = 'black'),
                 strip.text.y    = element_text(size = 17, colour = 'black'),
                 plot.caption    = element_text(size = 15, hjust = 0, colour = 'black'),
                 legend.position = 'bottom')
fig2
ggplot2::ggsave(filename = paste0('D:/Figure2_paper1.png'), plot = fig2, device = 'png', width = 6, height = 6.1, units = 'in', dpi = 350); rm(fig2)

## Figure 3
## Bivariate map of SPEI severity vs crop classes diversity (agriculture exposure)
# Crop classes diversity
crp_diversity <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_vulnerability/crop_classes_diversity.tif'))

## Index characterization
# Extreme trend
idx_sxt <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_vulnerability/SPEI-6_extreme_trend.tif'))
# # Extreme years' count
# idx_cnt <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_vulnerability/SPEI-6_extreme_years_count.tif'))

# Zonal statistics per extreme drought cluster (mean)
agex_sgn_extreme_trend <- terra::zonal(x = idx_sxt, z = agex_sgn_num, fun = 'mean', na.rm = T, as.raster = T)
agex_sgn_crp_diversity <- terra::zonal(x = crp_diversity, z = agex_sgn_num, fun = 'mean', na.rm = T, as.raster = T)

# Merge data.frames
agex_sgn_severity_vs_crops <- terra::as.data.frame(x = c(agex_sgn, agex_sgn_extreme_trend, agex_sgn_crp_diversity), xy = T, cell = T)
rm(agex_sgn_extreme_trend, agex_sgn_crp_diversity)

# Produce bivariate maps
severity_brks <- quantile(x = agex_sgn_metrics$`SPEI.6_slope_95th`, probs = seq(0,1,1/4))
diversity_brks <- quantile(x = agex_sgn_metrics$crop_classes_diversity, probs = seq(0,1,1/4))
agex_sgn_severity_vs_crops$Severity_class <- cut(x = agex_sgn_severity_vs_crops$`SPEI-6_slope_95th`, breaks = severity_brks) |> as.numeric() |> as.character()
agex_sgn_severity_vs_crops$Diversity_class <- cut(x = agex_sgn_severity_vs_crops$crop_classes_diversity, breaks = diversity_brks) |> as.numeric() |> as.character()
agex_sgn_severity_vs_crops$bi_class <- paste0(agex_sgn_severity_vs_crops$Severity_class,'-',agex_sgn_severity_vs_crops$Diversity_class)
agex_sgn_severity_vs_crops <- agex_sgn_severity_vs_crops[,c('cell','x','y','extreme_cluster','SPEI-6_slope_95th','crop_classes_diversity','bi_class')]

colours <- data.frame(bi_class = sort(unique(agex_sgn_severity_vs_crops$bi_class)))
colours$category <- 1:nrow(colours)

agex_sgn_severity_vs_crops <- dplyr::left_join(x = agex_sgn_severity_vs_crops, y = colours, by = 'bi_class')

aux <- agex_sgn
terra::values(aux) <- NA
aux[agex_sgn_severity_vs_crops$cell] <- agex_sgn_severity_vs_crops$category
names(aux) <- 'category'

make_legend <- function(dims = 4,
                        xlabl = 'Drought trend severity',
                        ylabl = 'Crop classes diversity',
                        xbrks = c('Negative','Low','Medium','High'),
                        ybrks = c('Very low','Low','Medium','High'),
                        palette = 'GrPink2',
                        categories = colours)
{
  categories$x <- strsplit(x = categories$bi_class, split = '-') |> purrr::map(1) |> unlist() |> as.numeric()
  categories$y <- strsplit(x = categories$bi_class, split = '-') |> purrr::map(2) |> unlist() |> as.numeric()
  legend <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = categories, mapping = ggplot2::aes(x = x, y = y, fill = bi_class)) +
    biscale::bi_scale_fill(pal = palette, dim = dims) +
    ggplot2::theme_minimal() +
    ggplot2::coord_equal() +
    ggplot2::xlab(xlabl) +
    ggplot2::ylab(ylabl) +
    ggplot2::scale_x_continuous(breaks = 1:dims, labels = xbrks) +
    ggplot2::scale_y_continuous(breaks = 1:dims, labels = ybrks) +
    ggplot2::theme(legend.position = 'none',
                   axis.text       = element_text(family = 'serif'),
                   axis.text.x     = element_text(size = 15, colour = 'black', angle = 45),
                   axis.text.y     = element_text(size = 15, colour = 'black', angle = 45),
                   axis.title      = element_text(size = 17, colour = 'black', face = 'bold', family = 'serif'),
                   axis.line       = element_blank(),
                   axis.ticks      = element_blank())
  return(legend)
}

ggm <- list()
for(i in 1:length(shp)){
  
  # Coordinates
  aux_dfm <- terra::as.data.frame(x = terra::crop(aux, terra::ext(shp[[i]])), xy = T, cell = T)
  aux_dfm <- dplyr::left_join(x = aux_dfm, y = colours, by = 'category')
  # Map
  ggm[[i]] <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = aux_dfm, aes(x, y, fill = bi_class)) +
    ggplot2::geom_sf(data = sf::st_as_sf(shp[[i]]), fill = NA) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = 'none', legend.key.width = unit(1.3, 'cm')) +
    biscale::bi_scale_fill(pal = 'GrPink2', dim = 4)
  
}

# fig3 <- ggpubr::ggarrange(ggm[[1]], ggm[[2]], ggm[[3]], ggm[[4]], ggm[[5]], ggm[[6]],
#                           labels = c('a)', 'b)', 'c)', 'd)', 'e)', 'f)'),
#                           ncol = 2, nrow = 3, font.label = list(size = 30, family = 'serif', face = 'plain'))
# 
# ggplot2::ggsave(filename = paste0('D:/Figure3_paper1.png'), plot = fig3, device = 'png', width = 10, height = 12.5, units = 'in', dpi = 350)

fig3_leg <- make_legend(dims = 4,
                        xlabl = 'Drought severity',
                        ylabl = 'Crop classes diversity',
                        xbrks = c('Negative','Low','Medium','High'),
                        ybrks = c('Very low','Low','Medium','High'),
                        palette = 'GrPink2',
                        categories = colours)
# ggplot2::ggsave(filename = paste0('D:/Figure3_paper1_legend.png'), plot = fig3_leg, device = 'png', width = 5, height = 5, units = 'in', dpi = 350)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Panel a) Cluster: 157 (USA)
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

idx_ts <- terra::rast('//CATALOGUE/WFP_ClimateRiskPr1/agroclimExtremes/agex_indices/agex_spei-6/agex_spei-6_25km/one_s1_spei-6_25km.tif')
# Filter time series per cluster of interest
aux <- agex_sgn; aux[aux != 157] <- NA
idx_ts_msk <- terra::mask(x = idx_ts, mask = aux) |> terra::trim()
# Time series plotting
idx_ts_msk_dfm <- terra::as.data.frame(x = idx_ts_msk, xy = T, cell = F)
gg_ts_NA <- idx_ts_msk_dfm |>
  tidyr::pivot_longer(cols = 3:ncol(idx_ts_msk_dfm), names_to = 'Years', values_to = 'SPEI-6') |>
  dplyr::mutate(Years = as.numeric(gsub('spei-6_', '', Years)),
                ID = paste0(x, '-', y)) |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = Years, y = `SPEI-6`, group = ID)) +
  ggplot2::geom_line(alpha = 0.1) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'bold', family = 'serif'))

# Add polygon extent to the corresponding map
ext_vct1 <- as.vector(terra::ext(idx_ts_msk))
ggm1 <- ggm[[1]] +
  ggplot2::geom_rect(aes(xmin = ext_vct1[1],
                         xmax = ext_vct1[2],
                         ymin = ext_vct1[3],
                         ymax = ext_vct1[4]),
                     color = 'red', size = 1.5, fill = NA)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Panel b) Cluster: 420 (Brazil, Paraguay) two seasons: S1
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

idx_ts <- terra::rast('//CATALOGUE/WFP_ClimateRiskPr1/agroclimExtremes/agex_indices/agex_spei-6/agex_spei-6_25km/two_s1_spei-6_25km.tif')
# Filter time series per cluster of interest
aux <- agex_sgn; aux[aux != 420] <- NA
idx_ts_msk <- terra::mask(x = idx_ts, mask = aux) |> terra::trim()
# Time series plotting
idx_ts_msk_dfm <- terra::as.data.frame(x = idx_ts_msk, xy = T, cell = F)
gg_ts_SA <- idx_ts_msk_dfm |>
  tidyr::pivot_longer(cols = 3:ncol(idx_ts_msk_dfm), names_to = 'Years', values_to = 'SPEI-6') |>
  dplyr::mutate(Years = as.numeric(gsub('spei-6_', '', Years)),
                ID = paste0(x, '-', y)) |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = Years, y = `SPEI-6`, group = ID)) +
  ggplot2::geom_line(alpha = 0.1) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'bold', family = 'serif'))

# Add polygon extent to the corresponding map
ext_vct2 <- as.vector(terra::ext(idx_ts_msk))
ggm2 <- ggm[[2]] +
  ggplot2::geom_rect(aes(xmin = ext_vct2[1],
                         xmax = ext_vct2[2],
                         ymin = ext_vct2[3],
                         ymax = ext_vct2[4]),
                     color = 'red', size = 1.5, fill = NA)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Panel c) Cluster: 270 (Ivory Coast, Mali, Guinea)
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

idx_ts <- terra::rast('//CATALOGUE/WFP_ClimateRiskPr1/agroclimExtremes/agex_indices/agex_spei-6/agex_spei-6_25km/one_s1_spei-6_25km.tif')
# Filter time series per cluster of interest
aux <- agex_sgn; aux[aux != 270] <- NA
idx_ts_msk <- terra::mask(x = idx_ts, mask = aux) |> terra::trim()
# Time series plotting
idx_ts_msk_dfm <- terra::as.data.frame(x = idx_ts_msk, xy = T, cell = F)
gg_ts_AF <- idx_ts_msk_dfm |>
  tidyr::pivot_longer(cols = 3:ncol(idx_ts_msk_dfm), names_to = 'Years', values_to = 'SPEI-6') |>
  dplyr::mutate(Years = as.numeric(gsub('spei-6_', '', Years)),
                ID = paste0(x, '-', y)) |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = Years, y = `SPEI-6`, group = ID)) +
  ggplot2::geom_line(alpha = 0.1) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'bold', family = 'serif'))

# Add polygon extent to the corresponding map
ext_vct3 <- as.vector(terra::ext(idx_ts_msk))
ggm3 <- ggm[[3]] +
  ggplot2::geom_rect(aes(xmin = ext_vct3[1],
                         xmax = ext_vct3[2],
                         ymin = ext_vct3[3],
                         ymax = ext_vct3[4]),
                     color = 'red', size = 1.5, fill = NA)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Panel d) Extreme Drought Cluster: 33 (Netherlands, others)
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

idx_ts <- terra::rast('//CATALOGUE/WFP_ClimateRiskPr1/agroclimExtremes/agex_indices/agex_spei-6/agex_spei-6_25km/one_s1_spei-6_25km.tif')
# Filter time series per cluster of interest
aux <- agex_sgn; aux[aux != 33] <- NA
idx_ts_msk <- terra::mask(x = idx_ts, mask = aux) |> terra::trim()
# Time series plotting
idx_ts_msk_dfm <- terra::as.data.frame(x = idx_ts_msk, xy = T, cell = F)
gg_ts_EU <- idx_ts_msk_dfm |>
  tidyr::pivot_longer(cols = 3:ncol(idx_ts_msk_dfm), names_to = 'Years', values_to = 'SPEI-6') |>
  dplyr::mutate(Years = as.numeric(gsub('spei-6_', '', Years)),
                ID = paste0(x, '-', y)) |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = Years, y = `SPEI-6`, group = ID)) +
  ggplot2::geom_line(alpha = 0.1) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'bold', family = 'serif'))

# Add polygon extent to the corresponding map
ext_vct4 <- as.vector(terra::ext(idx_ts_msk))
ggm4 <- ggm[[4]] +
  ggplot2::geom_rect(aes(xmin = ext_vct4[1],
                         xmax = ext_vct4[2],
                         ymin = ext_vct4[3],
                         ymax = ext_vct4[4]),
                     color = 'red', size = 1.5, fill = NA)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Panel e) Extreme Drought Cluster: 201 (India)
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

idx_ts <- terra::rast('//CATALOGUE/WFP_ClimateRiskPr1/agroclimExtremes/agex_indices/agex_spei-6/agex_spei-6_25km/two_s1_spei-6_25km.tif')
# Filter time series per cluster of interest
aux <- agex_sgn; aux[aux != 201] <- NA
idx_ts_msk <- terra::mask(x = idx_ts, mask = aux) |> terra::trim()
# Time series plotting
idx_ts_msk_dfm <- terra::as.data.frame(x = idx_ts_msk, xy = T, cell = F)
gg_ts_AS <- idx_ts_msk_dfm |>
  tidyr::pivot_longer(cols = 3:ncol(idx_ts_msk_dfm), names_to = 'Years', values_to = 'SPEI-6') |>
  dplyr::mutate(Years = as.numeric(gsub('spei-6_', '', Years)),
                ID = paste0(x, '-', y)) |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = Years, y = `SPEI-6`, group = ID)) +
  ggplot2::geom_line(alpha = 0.1) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'bold', family = 'serif'))

# Add polygon extent to the corresponding map
ext_vct5 <- as.vector(terra::ext(idx_ts_msk))
ggm5 <- ggm[[5]] +
  ggplot2::geom_rect(aes(xmin = ext_vct5[1],
                         xmax = ext_vct5[2],
                         ymin = ext_vct5[3],
                         ymax = ext_vct5[4]),
                     color = 'red', size = 1.5, fill = NA)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Panel f) Extreme Drought Cluster: 442 (Australia)
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

idx_ts <- terra::rast('//CATALOGUE/WFP_ClimateRiskPr1/agroclimExtremes/agex_indices/agex_spei-6/agex_spei-6_25km/one_s1_spei-6_25km.tif')
# Filter time series per cluster of interest
aux <- agex_sgn; aux[aux != 442] <- NA
idx_ts_msk <- terra::mask(x = idx_ts, mask = aux) |> terra::trim()
# Time series plotting
idx_ts_msk_dfm <- terra::as.data.frame(x = idx_ts_msk, xy = T, cell = F)
gg_ts_OC <- idx_ts_msk_dfm |>
  tidyr::pivot_longer(cols = 3:ncol(idx_ts_msk_dfm), names_to = 'Years', values_to = 'SPEI-6') |>
  dplyr::mutate(Years = as.numeric(gsub('spei-6_', '', Years)),
                ID = paste0(x, '-', y)) |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = Years, y = `SPEI-6`, group = ID)) +
  ggplot2::geom_line(alpha = 0.1) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
                 axis.text       = element_text(family = 'serif'),
                 axis.text.x     = element_text(size = 15, colour = 'black'),
                 axis.text.y     = element_text(size = 15, colour = 'black'),
                 axis.title      = element_text(size = 17, colour = 'black', face = 'bold', family = 'serif'))

# Add polygon extent to the corresponding map
ext_vct6 <- as.vector(terra::ext(idx_ts_msk))
ggm6 <- ggm[[6]] +
  ggplot2::geom_rect(aes(xmin = ext_vct6[1],
                         xmax = ext_vct6[2],
                         ymin = ext_vct6[3],
                         ymax = ext_vct6[4]),
                     color = 'red', size = 1.5, fill = NA)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Full figure
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

layout <- matrix(c(1:12,NA,13,13,NA), nrow = 4, ncol = 4, byrow = TRUE) # layout <- matrix(c(NA,1:2,NA,NA,3:4,5,NA,6:7,NA,NA,8,8,NA), nrow = 4, ncol = 4, byrow = TRUE)

fg_01 <- ggpubr::annotate_figure(gg_ts_NA, top = text_grob(label = 'a) Time series', face = 'plain', size = 15, family = 'serif'))
fg_02 <- ggpubr::annotate_figure(ggm1, top = text_grob(label = 'a)', face = 'plain', size = 20, family = 'serif'))
fg_03 <- ggpubr::annotate_figure(ggm2, top = text_grob(label = 'b)', face = 'plain', size = 20, family = 'serif'))
fg_04 <- ggpubr::annotate_figure(gg_ts_SA, top = text_grob(label = 'b) Time series', face = 'plain', size = 15, family = 'serif'))

fg_05 <- ggpubr::annotate_figure(gg_ts_AF, top = text_grob(label = 'c) Time series', face = 'plain', size = 15, family = 'serif'))
fg_06 <- ggpubr::annotate_figure(ggm3, top = text_grob(label = 'c)', face = 'plain', size = 20, family = 'serif'))
fg_07 <- ggpubr::annotate_figure(ggm4, top = text_grob(label = 'd)', face = 'plain', size = 20, family = 'serif'))
fg_08 <- ggpubr::annotate_figure(gg_ts_EU, top = text_grob(label = 'd) Time series', face = 'plain', size = 15, family = 'serif'))

fg_09 <- ggpubr::annotate_figure(gg_ts_AS, top = text_grob(label = 'e) Time series', face = 'plain', size = 15, family = 'serif'))
fg_10 <- ggpubr::annotate_figure(ggm5, top = text_grob(label = 'e)', face = 'plain', size = 20, family = 'serif'))
fg_11 <- ggpubr::annotate_figure(ggm6, top = text_grob(label = 'f)', face = 'plain', size = 20, family = 'serif'))
fg_12 <- ggpubr::annotate_figure(gg_ts_OC, top = text_grob(label = 'f) Time series', face = 'plain', size = 15, family = 'serif'))

fg_13 <- ggpubr::annotate_figure(fig3_leg, top = NULL)

all <- gridExtra::grid.arrange(fg_01, fg_02, fg_03, fg_04,
                               fg_05, fg_06, fg_07, fg_08,
                               fg_09, fg_10, fg_11, fg_12,
                               fg_13,
                               layout_matrix = layout)
ggplot2::ggsave(filename = paste0('D:/Figure3_paper1_test.png'), plot = all, device = 'png', width = 12, height = 10, units = 'in', dpi = 350)

### Figure 4
# Compute diversity of livestock
lvstck_units  <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_vulnerability/lsu_total.tif'))
agex_sgn_lvstck_units  <- terra::zonal(x = lvstck_units, z = agex_sgn_num, fun = 'mean', na.rm = T, as.raster = T)
agex_sgn_severity_vs_lvstck <- terra::as.data.frame(x = c(agex_sgn, agex_sgn_extreme_trend, agex_sgn_lvstck_units), xy = T, cell = T)

# Produce bivariate maps
diversity_brks <- quantile(x = agex_sgn_metrics$total_livestock_units, probs = seq(0,1,1/4))

lvs_dir <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory/data/_global/livestock'
anmls   <- list.dirs(path = lvs_dir, full.names = F, recursive = F)
anmls   <- anmls[-grep('horses',anmls)]
lvstc_fls <- list.files2(path = paste0(lvs_dir,'/',anmls), pattern = '_Da.tif$', full.names = T); rm(lvs_dir, anmls)
names(lvstc_fls) <- NULL
lvstc_cnt <- terra::rast(lvstc_fls)

# Computing Livestock Units
# LU: Livestock Unit
# LU = Buffaloes * 1 + Cattle * 1 + Chickens * ((0.007 + 0.014)/2) +
#      Ducks * 0.01 + Goats * 0.1 +Pigs * ((0.5+0.027)/2) + Sheep * 0.1
# Source: https://ec.europa.eu/eurostat/statistics-explained/index.php?title=Glossary:Livestock_unit_(LSU)
lsu <- lvstc_cnt[[1]] + lvstc_cnt[[2]] + lvstc_cnt[[3]]*((0.007 + 0.014)/2) + lvstc_cnt[[4]]*0.01 + lvstc_cnt[[5]]*0.1 + lvstc_cnt[[6]]*((0.5+0.027)/2) + lvstc_cnt[[7]]*0.1
rm(lvstc_fls, lvstc_cnt)
names(lsu) <- 'livestock_units'
lsu_25km <- terra::resample(x = lsu, y = agex_sgn, method = 'cubicspline', threads = T)
lsu_25km <- terra::mask(x = lsu_25km, mask = agex_sgn); rm(lsu)

sgn_lsu_diversity <- terra::zonal(x = lsu_25km, z = ref, fun = 'median', na.rm = T, as.raster = T)
names(sgn_lsu_diversity) <- 'livestock_diversity'

dfm <- terra::as.data.frame(x = c(agex_sgn,
                                  sgn_intensity, sgn_cnt_trend, sgn_ext_trend,
                                  sgn_lsu_diversity), xy = T, cell = T)

# Produce bivariate maps
dims <- 4
dfm_livestock_vs_intensity <- biscale::bi_class(.data = dfm, x = 'SPEI-6_extreme_trend', y = 'livestock_diversity', style = 'quantile', dim = dims) |> base::as.data.frame()
dfm_livestock_vs_intensity <- dfm_livestock_vs_intensity[,c('cell','x','y','extreme_cluster','SPEI-6_extreme_trend','livestock_diversity','bi_class')]

colours <- data.frame(bi_class = sort(unique(dfm_livestock_vs_intensity$bi_class)))
colours$category <- 1:nrow(colours)

dfm_livestock_vs_intensity <- dplyr::left_join(x = dfm_livestock_vs_intensity, y = colours, by = 'bi_class')

aux <- agex_sgn
terra::values(aux) <- NA
aux[dfm_livestock_vs_intensity$cell] <- dfm_livestock_vs_intensity$category
names(aux) <- 'category'

ggm <- list()
for(i in 1:length(shp)){
  
  # Coordinates
  aux_dfm <- terra::as.data.frame(x = terra::crop(aux, terra::ext(shp[[i]])), xy = T, cell = T)
  aux_dfm <- dplyr::left_join(x = aux_dfm, y = colours, by = 'category')
  # Map
  ggm[[i]] <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = aux_dfm, aes(x, y, fill = bi_class)) +
    ggplot2::geom_sf(data = sf::st_as_sf(shp[[i]]), fill = NA) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = 'none', legend.key.width = unit(1.3, 'cm')) +
    biscale::bi_scale_fill(pal = 'GrPink2', dim = dims)
  
}

fig3 <- ggpubr::ggarrange(ggm[[1]], ggm[[2]], ggm[[3]], ggm[[4]], ggm[[5]], ggm[[6]],
                          labels = c('a)', 'b)', 'c)', 'd)', 'e)', 'f)'),
                          ncol = 2, nrow = 3, font.label = list(size = 30, family = 'serif', face = 'plain'))
ggplot2::ggsave(filename = paste0('D:/Figure3_paper1.png'), plot = fig3, device = 'png', width = 10, height = 12.5, units = 'in', dpi = 350)



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

## Graph of global cluster One growing season
# World shapefile
wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
wrl <- terra::aggregate(wrl)
# Load extreme weather signatures
fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_results/clusters'), pattern = paste0('agex_global_',index,'_',gs,'_s',season,'_fmadogram_*.*.tif$'), full.names = T)
agex_sgn <- terra::rast(fls)
names(agex_sgn) <- 'extreme_signature'
# Colour palette
n_col <- length(unique(terra::values(agex_sgn,na.rm = T)))
col_pltt <- MetBrewer::met.brewer(name = 'Signac', n = n_col)
set.seed(1); col_pltt <- sample(x = col_pltt, size = n_col, replace = F); rm(n_col)
agex_sgn <- terra::trim(x = agex_sgn)
# Figure
png(filename = paste0('D:/agex_global_',index,'_',gs,'_s',season,'_fmadogram.png'), width = 3132, height = 2359, units = 'px')
plot(wrl, ext = terra::ext(agex_sgn))
plot(agex_sgn, add = T, col = col_pltt, plg = list(cex = 5), cex.main = 7)
dev.off()

## Time series figure
# Load extreme weather signatures
fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_results/clusters'), pattern = paste0('agex_global_',index,'_',gs,'_s',season,'_fmadogram_*.*.tif$'), full.names = T)
agex_sgn <- terra::rast(fls)
names(agex_sgn) <- 'extreme_signature'
# Load indices values
idx <- terra::rast(paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_25km/',gs,'_s',season,'_',index,'_25km.tif'))
idx_tsr <- terra::zonal(x = idx, z = agex_sgn, fun = 'median')

idx_tsr_lng <- idx_tsr |>
  tidyr::pivot_longer(cols = 2:ncol(idx_tsr), names_to = 'year', values_to = 'spei') |>
  base::as.data.frame()
idx_tsr_lng$year <- gsub('spei-6_','',idx_tsr_lng$year) |> as.numeric()
idx_tsr_lng$extreme_signature <- as.factor(idx_tsr_lng$extreme_signature)

tsr <- idx_tsr_lng |>
  ggplot2::ggplot(aes(x = year, y = spei, colour = extreme_signature)) +
  ggplot2::geom_line(alpha = 0.4) +
  ggplot2::scale_color_manual(values = col_pltt) +
  ggplot2::xlab('Year') +
  ggplot2::ylab('SPEI-6 median across signature') +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = 'blue') +
  ggplot2::theme_bw() +
  ggplot2::theme(text            = element_text(size = 17, colour = 'black'),
                 axis.text       = element_text(size = 16, colour = 'black'),
                 axis.title      = element_text(size = 20, colour = 'black'),
                 legend.text     = element_text(size = 13, colour = 'black'),
                 legend.title    = element_blank(),
                 plot.title      = element_text(size = 25, colour = 'black'),
                 plot.subtitle   = element_text(size = 17, colour = 'black'),
                 strip.text.x    = element_text(size = 17, colour = 'black'),
                 strip.text.y    = element_text(size = 17, colour = 'black'),
                 plot.caption    = element_text(size = 15, hjust = 0, colour = 'black'),
                 legend.position = 'none')
ggplot2::ggsave(filename = 'D:/time_series_per_signature.png', plot = tsr, device = 'png', units = 'in', width = 12, height = 8, dpi = 350)
