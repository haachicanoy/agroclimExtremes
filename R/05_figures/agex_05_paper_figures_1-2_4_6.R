## ------------------------------------------ ##
## Paper figures 1, 2, 4 & 6
## By: Harold Achicanoy
## WUR & ABC
## Oct. 2024
## ------------------------------------------ ##

## R options and packages loading ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,MetBrewer,rnaturalearth,tidyverse,ggpubr,tseries,biscale,pals,cowplot,gridExtra,grid))

# Activate Windows fonts
extrafont::font_import()
extrafont::loadfonts(device = 'win')
extrafont::fonts()

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
wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
wrl <- wrl[,c('adm0_iso','name_en','continent','subregion')]

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

## Figure 1 ----
# Extreme drought clusters per continent
# Color palette definition
n_col    <- length(unique(terra::values(agex_sgn,na.rm = T))) # Number of unique colors needed
col_pltt <- MetBrewer::met.brewer(name = 'Signac', n = n_col) # Color palette
set.seed(1235); col_pltt <- sample(x = col_pltt, size = n_col, replace = F); rm(n_col) # Randomize colors assignment

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
  
}; rm(i,stp,sub_stp,n_col,col_pltt)

# Organize maps into one graph adding panel labels (alternative 1)
# fig1 <- ggpubr::ggarrange(ggm[[1]], ggm[[2]], ggm[[3]], ggm[[4]], ggm[[5]], ggm[[6]],
#                              labels = c('a)', 'b)', 'c)', 'd)', 'e)', 'f)'),
#                              ncol = 2, nrow = 3, font.label = list(size = 30, family = 'serif', face = 'plain')); fig1

# Organize maps into one graph adding panel labels (alternative 2)
fg_01 <- ggpubr::annotate_figure(ggm[[1]], top = text_grob(label = 'a) North America', face = 'plain', size = 20, family = 'serif'))
fg_02 <- ggpubr::annotate_figure(ggm[[2]], top = text_grob(label = 'b) South America', face = 'plain', size = 20, family = 'serif'))
fg_03 <- ggpubr::annotate_figure(ggm[[3]], top = text_grob(label = 'c) Africa', face = 'plain', size = 20, family = 'serif'))
fg_04 <- ggpubr::annotate_figure(ggm[[4]], top = text_grob(label = 'd) Europe', face = 'plain', size = 20, family = 'serif'))
fg_05 <- ggpubr::annotate_figure(ggm[[5]], top = text_grob(label = 'e) Asia', face = 'plain', size = 20, family = 'serif'))
fg_06 <- ggpubr::annotate_figure(ggm[[6]], top = text_grob(label = 'f) Oceania', face = 'plain', size = 20, family = 'serif'))

layout <- matrix(1:6, nrow = 3, ncol = 2, byrow = T)

fig1 <- gridExtra::grid.arrange(fg_01, fg_02,
                                fg_03, fg_04,
                                fg_05, fg_06,
                                layout_matrix = layout); rm(ggm,fg_01,fg_02,fg_03,fg_04,fg_05,fg_06,layout)
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/Figure1_paper1.png'), plot = fig1, device = 'png', width = 10, height = 12.5, units = 'in', dpi = 350); rm(fig1)

## Figure 2 ----
# Percentage of extreme drought clusters in one, two, or more countries (Collaboration potential)
# Color palette definition
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
  ggplot2::xlab('Number of countries per extreme drought cluster') +
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
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/Figure2_paper1.png'), plot = fig2, device = 'png', width = 7, height = 6.1, units = 'in', dpi = 350); rm(fig2,col_pltt)

## Figure 4 ----
crops_hotspots <- utils::read.csv(paste0(root,'/agroclimExtremes/agex_results/agex_crops_hotspots.csv'))
crops_hotspots$continents[crops_hotspots$continents == 'Asia,Europe'] <- 'Asia'

## Figure 4
fig4 <- crops_hotspots |>
  ggplot2::ggplot(aes(x = reorder(continents, SPEI.6_slope_95th, FUN = median, decreasing = T), y = SPEI.6_slope_95th * 10)) +
  ggplot2::geom_boxplot() +
  ggplot2::ylim(0, 0.6) +
  ggplot2::xlab('Continent') +
  ggplot2::ylab('SPEI units per decade') +
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
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/Figure4_paper1.png'), plot = fig4, device = 'png', width = 7, height = 6.1, units = 'in', dpi = 350); rm(fig4)

## Figure 6 ----
mtx_ref <- expand.grid(1:4, 1:4)
names(mtx_ref) <- c('Severity','Diversity')
mtx_ref <- mtx_ref |>
  dplyr::arrange(Severity) |>
  dplyr::mutate(key = paste0(Severity,'-',Diversity)) |>
  base::as.data.frame()
livestock_hotspots <- utils::read.csv(paste0(root,'/agex_results/agex_all_metrics_hotspots.csv'))
livestock_hotspots <- livestock_hotspots |>
  dplyr::left_join(y = mtx_ref,
                   by = c('livestock_diversity_exposure' = 'key')) |>
  dplyr::select(extreme_cluster,growing_seasons,
                countries_count,countries,continents,
                Severity,Diversity,SPEI.6_slope_95th) |>
  dplyr::filter(Severity >= 3 & Diversity >= 3)
livestock_hotspots$continents[grep(pattern = ',', x = livestock_hotspots$continents)] <- c('Asia','Europe','Europe','Asia','Europe','Africa','Asia','Asia','Asia','Asia','Asia','Asia','Europe')
livestock_hotspots$continents <- factor(x = livestock_hotspots$continents, levels = c('Africa','Asia','South America','Europe','North America'))

fig6 <- livestock_hotspots |>
  ggplot2::ggplot(aes(x = continents, y = SPEI.6_slope_95th * 10)) +
  ggplot2::geom_boxplot() +
  ggplot2::ylim(0, 0.5) +
  ggplot2::xlab('Continent') +
  ggplot2::ylab('SPEI units per decade') +
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
ggplot2::ggsave(filename = paste0(root,'/agex_results/agex_figures/Figure6_paper1.png'), plot = fig6, device = 'png', width = 7, height = 6.1, units = 'in', dpi = 350); rm(fig6)
