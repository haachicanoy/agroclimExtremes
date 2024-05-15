## ------------------------------------------ ##
## Paper figures/maps
## By: Harold Achicanoy
## WUR & ABC
## May 2024
## ------------------------------------------ ##

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,MetBrewer,rnaturalearth,tidyverse,ggpubr,tseries,biscale,pals,cowplot,gridExtra))

list.files2 <- Vectorize(FUN = list.files, vectorize.args = 'path')
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')
get_trend <- function(x){
  if(!all(is.na(x))){
    x <- as.numeric(na.omit(x))
    y <- as.numeric(trend::sens.slope(x)$estimates)
  } else { y <- NA }
  return(y)
}
get_ext_trend <- function(x){
  if(!all(is.na(x))){
    x <- as.numeric(na.omit(x))
    dfm <- data.frame(time = 1:length(x), ts = x)
    qrf <- quantreg::rq(formula = ts ~ time, tau = .90, data = dfm)
    y <- as.numeric(qrf$coefficients[2])
  } else { y <- NA }
  return(y)
}
getmode <- function(v){
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
cv_fun <- function(x, na.rm = T){
  sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm) * 100
}

## Key arguments
root   <- '//CATALOGUE/WFP_ClimateRiskPr1'
index  <- 'spei-6'
# gs     <- 'one'
# season <- 1

### Extreme clusters (categorical)
agex_sgn <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn) <- 'extreme_cluster'
cls <- data.frame(extreme_cluster = sort(unique(terra::values(agex_sgn))))
cls$id <- cls$extreme_cluster
cls <- cls[,c('id','extreme_cluster')]
cls$extreme_cluster <- as.character(cls$extreme_cluster)
levels(agex_sgn) <- cls

### Load global shapefile
wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
wrl <- wrl[,c('adm0_iso','name_en','continent','subregion')]

afr <- rnaturalearth::ne_countries(scale = 'large', continent = 'africa', returnclass = 'sv')
eur <- rnaturalearth::ne_countries(scale = 'large', continent = 'europe', returnclass = 'sv')
eur <- terra::crop(x = eur, y = terra::ext(c(-20,100,30,75))) # xmin, xmax = 150, ymin, ymax
asi <- rnaturalearth::ne_countries(scale = 'large', continent = 'asia', returnclass = 'sv')
oce <- rnaturalearth::ne_countries(scale = 'large', continent = 'oceania', returnclass = 'sv')
oce <- terra::crop(x = oce, y = terra::ext(c(90,180,-54.75,20.55)))
nam <- rnaturalearth::ne_countries(scale = 'large', continent = 'north america', returnclass = 'sv')
nam <- terra::crop(x = nam, y = terra::ext(c(-140,-50,-0.39,60)))
sam <- rnaturalearth::ne_countries(scale = 'large', continent = 'south america', returnclass = 'sv')

shp <- list(nam, sam, afr, eur, asi, oce); rm(nam, sam, afr, eur, asi, oce)

# agex_sgn_afr <- terra::crop(x = agex_sgn, y = terra::ext(afr)); # agex_sgn_afr <- terra::mask(x = agex_sgn_afr, mask = afr)
# agex_sgn_eur <- terra::crop(x = agex_sgn, y = terra::ext(eur)); # agex_sgn_eur <- terra::mask(x = agex_sgn_eur, mask = eur)
# agex_sgn_asi <- terra::crop(x = agex_sgn, y = terra::ext(asi)); # agex_sgn_asi <- terra::mask(x = agex_sgn_asi, mask = asi)
# agex_sgn_oce <- terra::crop(x = agex_sgn, y = terra::ext(oce)); # agex_sgn_oce <- terra::mask(x = agex_sgn_oce, mask = oce)
# agex_sgn_nam <- terra::crop(x = agex_sgn, y = terra::ext(nam)); # agex_sgn_nam <- terra::mask(x = agex_sgn_nam, mask = nam)
# agex_sgn_sam <- terra::crop(x = agex_sgn, y = terra::ext(sam)); # agex_sgn_sam <- terra::mask(x = agex_sgn_sam, mask = sam)

# Color palette
n_col <- length(unique(terra::values(agex_sgn,na.rm = T)))
col_pltt <- MetBrewer::met.brewer(name = 'Signac', n = n_col)
set.seed(1); col_pltt <- sample(x = col_pltt, size = n_col, replace = F); rm(n_col)

## Figure 1
# myTheme <- rasterVis::rasterTheme(region = col_pltt)
# # Europe
# rasterVis::levelplot(x            = agex_sgn_eur,
#                      layers       = 'extreme_cluster',
#                      par.settings = myTheme,
#                      colorkey     = F,
#                      maxpixels    = terra::ncell(agex_sgn))

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

fig1 <- ggpubr::ggarrange(ggm[[1]], ggm[[2]], ggm[[3]], ggm[[4]], ggm[[5]], ggm[[6]],
                             labels = c('a)', 'b)', 'c)', 'd)', 'e)', 'f)'),
                             ncol = 2, nrow = 3, font.label = list(size = 30, family = 'serif', face = 'plain'))
ggplot2::ggsave(filename = paste0('D:/Figure1_paper1.png'), plot = fig1, device = 'png', width = 10, height = 12.5, units = 'in', dpi = 350)

### Figure 2
## Crop types diversity (agriculture exposure)
crp_diversity <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_vulnerability/crop_classes_diversity.tif'))
lvstck_units  <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_vulnerability/lsu_total.tif'))

## Index characterization
# Extreme trend
idx_sxt <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_vulnerability/SPEI-6_extreme_trend.tif'))
# Extreme years' count
idx_cnt <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_vulnerability/SPEI-6_extreme_years_count.tif'))

# Extreme clusters (numerical)
ref <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(ref) <- 'extreme_cluster'

agex_sgn_extreme_trend <- terra::zonal(x = idx_sxt, z = ref, fun = 'mean', na.rm = T, as.raster = T)
agex_sgn_crp_diversity <- terra::zonal(x = crp_diversity, z = ref, fun = 'mean', na.rm = T, as.raster = T)
agex_sgn_lvstck_units  <- terra::zonal(x = lvstck_units, z = ref, fun = 'mean', na.rm = T, as.raster = T)

agex_sgn_severity_vs_crops  <- terra::as.data.frame(x = c(agex_sgn, agex_sgn_extreme_trend, agex_sgn_crp_diversity), xy = T, cell = T)
agex_sgn_severity_vs_lvstck <- terra::as.data.frame(x = c(agex_sgn, agex_sgn_extreme_trend, agex_sgn_lvstck_units), xy = T, cell = T)

dfm <- terra::zonal(x = c(idx_sxt, crp_diversity, lvstck_units), z = ref, fun = 'mean', na.rm = T)

# Produce bivariate maps
severity_brks <- quantile(x = dfm$`SPEI-6_slope_95th`, probs = seq(0,1,1/4))
diversity_brks <- quantile(x = dfm$crop_classes_diversity, probs = seq(0,1,1/4))
agex_sgn_severity_vs_crops$`SPEI-6_slope_95th__class` <- cut(x = agex_sgn_severity_vs_crops$`SPEI-6_slope_95th`, breaks = severity_brks) |> as.numeric() |> as.character()
agex_sgn_severity_vs_crops$crop_classes_diversity__class <- cut(x = agex_sgn_severity_vs_crops$crop_classes_diversity, breaks = diversity_brks) |> as.numeric() |> as.character()
agex_sgn_severity_vs_crops$bi_class <- paste0(agex_sgn_severity_vs_crops$`SPEI-6_slope_95th__class`,'-',agex_sgn_severity_vs_crops$crop_classes_diversity__class)
agex_sgn_severity_vs_crops <- agex_sgn_severity_vs_crops[,c('cell','x','y','extreme_cluster','SPEI-6_slope_95th','crop_classes_diversity','bi_class')]

# map <- ggplot2::ggplot() +
#   ggplot2::geom_tile(data = dfm_crops_vs_intensity, aes(x, y, fill = bi_class)) +
#   ggplot2::geom_sf(data = sf::st_as_sf(wrl), fill = NA) +
#   ggplot2::coord_sf() +
#   ggplot2::theme_void() +
#   ggplot2::theme(legend.position = 'none', legend.key.width = unit(1.3, 'cm')) +
#   biscale::bi_scale_fill(pal = 'GrPink2', dim = dims)

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

# legend <- bi_legend(pal = 'GrPink2',
#                     dim = dims,
#                     xlab = 'Higher drought intensity',
#                     ylab = 'Higher crop types diversity',
#                     size = 10)
# 
# finalPlot <- cowplot::ggdraw() + cowplot::draw_plot(map, 0, 0, 1, 1) + cowplot::draw_plot(legend, -.03, 0.12, 0.3, 0.3)
# finalPlot

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

fig2 <- ggpubr::ggarrange(ggm[[1]], ggm[[2]], ggm[[3]], ggm[[4]], ggm[[5]], ggm[[6]],
                          labels = c('a)', 'b)', 'c)', 'd)', 'e)', 'f)'),
                          ncol = 2, nrow = 3, font.label = list(size = 30, family = 'serif', face = 'plain'))

ggplot2::ggsave(filename = paste0('D:/Figure2_paper1.png'), plot = fig2, device = 'png', width = 10, height = 12.5, units = 'in', dpi = 350)
fig2_leg <- make_legend(dims = 4,
                        xlabl = 'Drought severity',
                        ylabl = 'Crop classes diversity',
                        xbrks = c('Negative','Low','Medium','High'),
                        ybrks = c('Very low','Low','Medium','High'),
                        palette = 'GrPink2',
                        categories = colours)
ggplot2::ggsave(filename = paste0('D:/Figure2_paper1_legend.png'), plot = fig2_leg, device = 'png', width = 5, height = 5, units = 'in', dpi = 350)

# Extreme cluster 33. Netherlands
idx_ts <- terra::rast('//CATALOGUE/WFP_ClimateRiskPr1/agroclimExtremes/agex_indices/agex_spei-6/agex_spei-6_25km/one_s1_spei-6_25km.tif')
idx_ts <- -1 * idx_ts

aux <- agex_sgn
aux[aux != 33] <- NA
idx_ts_msk <- terra::mask(x = idx_ts, mask = aux)
idx_ts_msk <- terra::trim(x = idx_ts_msk)

idx_ts_msk_dfm <- terra::as.data.frame(x = idx_ts_msk, xy = T, cell = F)
plot(x = 1980:2022, y = as.numeric(idx_ts_msk_dfm[1,3:ncol(idx_ts_msk_dfm)]), ty = 'l', xlab = 'Years', ylab = 'SPEI-6', ylim = c(-3,2))
for(i in 2:nrow(idx_ts_msk_dfm)){
  lines(x = 1980:2022, y = as.numeric(idx_ts_msk_dfm[i,3:ncol(idx_ts_msk_dfm)]), col = 'red')
}

ggts <- idx_ts_msk_dfm |>
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
  # ggplot2::geom_smooth(alpha = 0.9, se = F)
ggplot2::ggsave(filename = paste0('D:/Figure2_paper1_EU_ts.png'), plot = ggts, device = 'png', width = 6, height = 5, units = 'in', dpi = 350)

### Figure 3
# Compute diversity of livestock
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
