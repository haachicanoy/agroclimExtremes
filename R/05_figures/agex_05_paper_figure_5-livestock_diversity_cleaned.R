## ------------------------------------------ ##
## Paper figure 5: LSU diversity vs drought intensification
## By: Harold Achicanoy
## WUR & ABC
## Sep 2025
## ------------------------------------------ ##

## R options and packages loading ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,MetBrewer,sf,rnaturalearth,tidyverse,ggpubr,tseries,biscale,pals,cowplot,gridExtra,grid))

# Activate Windows fonts
extrafont::font_import()
extrafont::loadfonts(device = 'win')
extrafont::fonts()

# Functions
make_maps_legend <- function(dims = 4,
                             xlabl = 'Drought trend severity',
                             ylabl = 'Livestock species diversity',
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
                   axis.text.x     = element_text(size = 14, colour = 'black', angle = 45, hjust = 1),
                   axis.text.y     = element_text(size = 14, colour = 'black', angle = 45),
                   axis.title      = element_text(size = 16, colour = 'black', face = 'bold', family = 'serif'),
                   axis.line       = element_blank(),
                   axis.ticks      = element_blank())
  return(legend)
}

## Setup arguments ----
root   <- '//CATALOGUE/AgroclimExtremes'
index  <- 'spei-6'

## Extreme drought clusters ----
# Categorical
agex_sgn <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn) <- 'extreme_cluster'
cls <- data.frame(extreme_cluster = sort(unique(terra::values(agex_sgn))))
cls$id <- cls$extreme_cluster
cls <- cls[,c('id','extreme_cluster')]
cls$extreme_cluster <- as.character(cls$extreme_cluster)
levels(agex_sgn) <- cls; rm(cls)

# Numerical
agex_sgn_num <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn_num) <- 'extreme_cluster'

# Metrics
agex_sgn_metrics <- utils::read.csv(paste0(root,'/agex_results/agex_all_metrics.csv'))

## Shapefiles ----
# Continent shapefiles. Some of them cropped
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

## Figure 5-supplementary ----
# Bivariate map of SPEI severity vs livestock units diversity (livestock exposure)

# Crop classes diversity
lsu_diversity <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/lsu_diversity_25km.tif'))

# Index characterization: extreme trend
idx_sxt <- terra::rast(paste0(root,'/agex_results/agex_vulnerability/SPEI-6_extreme_trend.tif'))

# Zonal statistics per extreme drought cluster (mean)
agex_sgn_extreme_trend <- terra::zonal(x = idx_sxt, z = agex_sgn_num, fun = 'mean', na.rm = T, as.raster = T)
metrics1 <- terra::zonal(x = idx_sxt, z = agex_sgn_num, fun = 'mean', na.rm = T)
agex_sgn_lsu_diversity <- terra::zonal(x = lsu_diversity, z = agex_sgn_num, fun = 'mean', na.rm = T, as.raster = T)
metrics2 <- terra::zonal(x = lsu_diversity, z = agex_sgn_num, fun = 'mean', na.rm = T)
metrics <- dplyr::left_join(x = metrics1, y = metrics2, by = 'extreme_cluster'); rm(metrics1, metrics2)

# Merge data.frames
agex_sgn_severity_vs_lvst <- terra::as.data.frame(x = c(agex_sgn, agex_sgn_extreme_trend, agex_sgn_lsu_diversity), xy = T, cell = T)
rm(agex_sgn_extreme_trend, agex_sgn_lsu_diversity)

# Produce bivariate maps
# Define severity and diversity quantiles
severity_brks <- quantile(x = metrics$`SPEI-6_slope_95th`, probs = seq(0,1,1/4))
severity_brks[1] <- terra::as.data.frame(idx_sxt) |> min()
severity_brks[5] <- terra::as.data.frame(idx_sxt) |> max()
diversity_brks <- quantile(x = metrics$livestock_units_diversity, probs = seq(0,1,1/4))
diversity_brks[1] <- terra::as.data.frame(lsu_diversity)|> min()
diversity_brks[5] <- terra::as.data.frame(lsu_diversity)|> max()
# Create categories
agex_sgn_severity_vs_lvst$Severity_class <- cut(x = agex_sgn_severity_vs_lvst$`SPEI-6_slope_95th`, breaks = severity_brks) |> as.numeric() |> as.character()
agex_sgn_severity_vs_lvst$Diversity_class <- cut(x = agex_sgn_severity_vs_lvst$livestock_units_diversity, breaks = diversity_brks) |> as.numeric() |> as.character()
# Create bivariate categories
agex_sgn_severity_vs_lvst$bi_class <- paste0(agex_sgn_severity_vs_lvst$Severity_class,'-',agex_sgn_severity_vs_lvst$Diversity_class)
agex_sgn_severity_vs_lvst <- agex_sgn_severity_vs_lvst[,c('cell','x','y','extreme_cluster','SPEI-6_slope_95th','livestock_units_diversity','bi_class')]
# Define color palette for bivariate categories
colours <- data.frame(bi_class = sort(unique(agex_sgn_severity_vs_lvst$bi_class)))
colours$category <- 1:nrow(colours)
agex_sgn_severity_vs_lvst <- dplyr::left_join(x = agex_sgn_severity_vs_lvst, y = colours, by = 'bi_class')

# Create auxiliary raster of bivariate categories
aux <- agex_sgn
terra::values(aux) <- NA
aux[agex_sgn_severity_vs_lvst$cell] <- agex_sgn_severity_vs_lvst$category
names(aux) <- 'category'

# Create a list of maps (per continent)
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
  
}; rm(i,aux_dfm)

### Legend ----
fig7_leg <- make_maps_legend(dims = 4,
                             xlabl = 'Extreme drought\nintensification',
                             ylabl = 'Livestock species\ndiversity',
                             xbrks = c('Negative','Low','Medium','High'),
                             ybrks = c('Very low','Low','Medium','High'),
                             palette = 'GrPink2',
                             categories = colours) # ggplot2::ggsave(filename = paste0('D:/Figure3_paper1_legend.png'), plot = fig3_leg, device = 'png', width = 5, height = 5, units = 'in', dpi = 350)

ggm_NA <- ggm[[1]]
ggm_SA <- ggm[[2]]
ggm_AF <- ggm[[3]]
ggm_EU <- ggm[[4]]
ggm_AS <- ggm[[5]]
ggm_OC <- ggm[[6]]

## Full figure ----
layout <- matrix(c(1:6,7,7), nrow = 4, ncol = 2, byrow = T)

fg_01 <- ggpubr::annotate_figure(ggm_NA,
                                 top = text_grob(
                                   label  = 'a) North America',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_02 <- ggpubr::annotate_figure(ggm_SA,
                                 top = text_grob(
                                   label  = 'b) South America',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_03 <- ggpubr::annotate_figure(ggm_AF,
                                 top = text_grob(
                                   label  = 'c) Africa',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_04 <- ggpubr::annotate_figure(ggm_EU,
                                 top = text_grob(
                                   label  = 'd) Europe',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_05 <- ggpubr::annotate_figure(ggm_AS,
                                 top = text_grob(
                                   label  = 'e) Asia',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_06 <- ggpubr::annotate_figure(ggm_OC,
                                 top = text_grob(
                                   label  = 'f) Oceania',
                                   face   = 'plain',
                                   size   = 17,
                                   family = 'serif'
                                 ))
fg_07 <- ggpubr::annotate_figure(fig7_leg, top = NULL)

fig3 <- gridExtra::grid.arrange(fg_01, fg_02,
                                fg_03, fg_04,
                                fg_05, fg_06,
                                fg_07,
                                layout_matrix = layout)
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/Figure5_paper1-livestock_diversity_cleaned.png'), plot = fig3, device = 'png', width = 10, height = 12.5, units = 'in', dpi = 350)
rm(ggm,
   ggm_NA,ggm_SA,ggm_AF,ggm_EU,ggm_AS,ggm_OC,
   gg_ts_NA,gg_ts_SA,gg_ts_AF,gg_ts_EU,gg_ts_AS,gg_ts_OC,
   fg_01,fg_02,fg_03,fg_04,fg_05,fg_06,fg_07,fg_08,fg_09,fg_10,fg_11,fg_12,fg_13,
   fig3_leg,fig3)
